#' LONGO analysis
#'
#' This is a internal function. Used by the package to take a formatted
#' dataframe and analyzes it based on the user inputted values.
#'
#' @param analyzeData dataframe to be analyzed
#' @param multiProbes Option to select what to do with multiple probes to a
#' single gene, default=mean
#' @param windowSize Option to alter the size of the windows to be created
#' from the dataframe, default=200
#' @param stepSize Option to alter the size of the steps used to create the
#'  windows, default=40
#' @param windowStyle Option to alter the method used to create windows,
#' default=mean
#' @param filterData Option to filter data using edgeR::cpm, default=TRUE
#' @param normalizeData Option to normalize data using
#' preprocessCore::normalize.quantiles, default=TRUE
#' @param controlColumn column of the data used as a control for the
#' statistical analyses
#' @import edgeR
#' @import preprocessCore
#' @import data.table
#' @return returns the simplified data
analyze <- function(analyzeData,
                    multiProbes,
                    windowSize,
                    stepSize,
                    windowStyle,
                    filterData,
                    normalizeData,
                    controlColumn) {

    KLdist <- function(freqs1, freqs2, unit=c("log", "log2", "log10")) {
        unit <- match.arg(unit)
        freqs1 <- freqs1 / sum(freqs1)
        freqs2 <- freqs2 / sum(freqs2)

        LR <- ifelse(freqs1 > 0, log(freqs1 / freqs2), 0)
        KL <- sum(freqs1 * LR)
        if (unit == "log2")
            KL <- KL / log(2)  # change from log to log2 scale
        if (unit == "log10")
            KL <- KL / log(10) # change from log to log10 scale
        return(KL)
    }

    JSdist <- function(x, y){
        0.5 * KLdist(x, (x + y) / 2) + 0.5 * KLdist(y, (x + y) / 2)
    }

    analyzeData <- as.data.frame(analyzeData)
    for (x in 2:(ncol(analyzeData) - 2)) {
        analyzeData[, x] <- as.numeric(as.character(analyzeData[, x]))
    }

    # remove duplicate gene reads
    if (multiProbes == "highest") {
        temp1.df <- analyzeData
        temp1.df$rowsum <- rowSums(temp1.df[, 2:(ncol(temp1.df) - 2)])
        temp1.df <- temp1.df[order(temp1.df$rowsum, decreasing=TRUE,
                                    na.last=NA),]
        row.names(temp1.df) <- 1:nrow(temp1.df)
        temp1.df <- (temp1.df[!duplicated(temp1.df[, (ncol(temp1.df) - 2)]),])
        temp1.df$rowsum <- NULL
        analyzeData <- temp1.df[, c((ncol(temp1.df) - 1),
            2:(ncol(temp1.df) - 2), (ncol(temp1.df)))]
    }
    else if (multiProbes == "mean") {
        analyzeData <- data.table::as.data.table(analyzeData)
        temp2 <- analyzeData[, lapply(.SD, mean, na.rm= TRUE),
            by=c(colnames(analyzeData)[(ncol(analyzeData)-1)],
            colnames(analyzeData)[ncol(analyzeData)]),
            .SDcols=c(2:(ncol(analyzeData)-2))]
        temp2 <-  as.data.frame(temp2)
        analyzeData <- temp2[c(1,3:ncol(temp2),2)]
    }

    # remove genes whose length are less than the median
    geneMedian <- median(analyzeData[,ncol(analyzeData)])
    analyzeData <- analyzeData[(analyzeData[,ncol(analyzeData)] > geneMedian),]

    # format of data after above
    #  1st col #####  many cols ####    last col ###
    #  symbol    |     expression      | length

    if (filterData) {
        isexpr <- rowSums(edgeR::cpm(
            analyzeData[2:(ncol(analyzeData)-1)]) > 1) >= 4
        analyzeData <- analyzeData[isexpr,]
    }

    if (normalizeData) {
        analyzeData[2:(ncol(analyzeData) - 1)] <-
            preprocessCore::normalize.quantiles(as.matrix(
                analyzeData[2:(ncol(analyzeData) - 1)]))
    }

    analyzeData$length <- as.numeric(as.character(analyzeData[,ncol(analyzeData)]))
    analyzeData <-
        analyzeData[order(analyzeData$length, decreasing=FALSE, na.last=NA),]
    row.names(analyzeData) <- 1:nrow(analyzeData)


    if(nrow(analyzeData)<200){
        return(NULL)
    }
    # calculate number of steps needed for the data
    s <- ((nrow(analyzeData) - windowSize) / stepSize) + 1
    #create dataframe for results
    simp.df <- data.frame("kb_length"=rep(0, s))
    # get labels of expression reads
    labels <- c("kb_length", colnames(analyzeData)[2:(ncol(analyzeData) - 1)])
    simp.df[labels] <- 0
    # create dataframe for pvalues
    pValues.df <- data.frame("kb_length"=rep(0, s))
    pValues.df[labels] <- 0
    # create dataframe for kL values
    JSdistances.df <- data.frame("kb_length"=rep(0, s))
    JSdistances.df[labels] <- 0

    w <- 1
    stepNumber <- 1
    while (stepNumber < s) {
        # get average length
        count <- 0
        if (windowStyle == "mean") {
            count <- sum(analyzeData$length[(w:(w + windowSize - 1))])
            count <- count / windowSize
        }
        if (windowStyle == "median") {
            count <- median(analyzeData$length[(w:(w + windowSize - 1))])
        }
        simp.df[stepNumber, 1] <- count
        JSdistances.df[stepNumber, 1] <- count
        pValues.df[stepNumber, 1] <-  count
        for (columnIndex in (2:(ncol(analyzeData) - 1))) {
            count <- 0
            if (windowStyle == "mean") {
                count <- sum(analyzeData[w:(w + windowSize - 1), columnIndex])
                count <- count / windowSize
            }
            if (windowStyle == "median") {
                count <- median(analyzeData[w:(w + windowSize - 1),
                                        columnIndex])
            }
            simp.df[stepNumber, columnIndex] <- count

            x <- analyzeData[w:(w + windowSize - 1), controlColumn]
            y <- analyzeData[w:(w + windowSize - 1), columnIndex]
            pvalue <- stats::wilcox.test(x, y, exact=FALSE)$p.value
            if (!is.nan(pvalue)) {
                pValues.df[stepNumber, columnIndex] <- pvalue
            } else {
                pValues.df[stepNumber, columnIndex] <- 1
            }
        }
        w <- w + stepSize
        stepNumber <- stepNumber + 1
    }
    for(stepNumber in (2:s)){
        for (columnIndex in (2:(ncol(analyzeData) - 1))) {
            x <- simp.df[2:stepNumber, controlColumn]
            y <- simp.df[2:stepNumber, columnIndex]
            if (all(x == 0) && all(y == 0)) {
                JSdistance <- 0
            } else {
                JSdistance <- JSdist(x, y)
            }
            if (is.na(JSdistance)) {
                JSdistance <- 0
            }
            JSdistances.df[stepNumber, columnIndex] <- JSdistance
        }
    }

    LQ.df <- LONGOquotient(simp.df, JSdistances.df, controlColumn)
    return(list(simp.df, pValues.df, JSdistances.df, LQ.df))
}
