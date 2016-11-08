#' LONGO analysis
#'
#' This is a internal function. Used by the package to take a formatted dataframe and analyzes it based on the user inputted values.
#'
#' @param DATA.DF dataframe to be analyzed
#' @param MULTI_PROBES Option to select what to do with multiple probes to a single gene, default=highest
#' @param WINDOW_SIZE Option to alter the size of the windows to be created from the dataframe, default=200
#' @param STEP_SIZE Option to alter the size of the steps used to create the windows, default=40
#' @param WINDOW_STYLE Option to alter the method used to create windows, default=mean
#' @param FILTERED Option to filter data, default=TRUE
#' @param NORMALIZED Option to normalize data, default=TRUE
#' @param CONTROL_COLUMN_INDEX column of the data used as a control for the statistical analyses
#' @import edgeR
#' @import preprocessCore
#' @import data.table
#' @return returns the simplified data
analyze <-
  function(DATA.DF,
           MULTI_PROBES,
           WINDOW_SIZE,
           STEP_SIZE,
           WINDOW_STYLE,
           FILTERED,
           NORMALIZED,
           CONTROL_COLUMN_INDEX) {
    KLdist <-
      function(freqs1,
               freqs2,
               unit = c("log", "log2", "log10"))
      {
        unit = match.arg(unit)
        freqs1 = freqs1 / sum(freqs1) # just to make sure ...
        freqs2 = freqs2 / sum(freqs2) # just to make sure ...

        if (any(!(freqs2 > 0)))
          warning("Vanishing value(s) in argument freqs2!")

        LR = ifelse(freqs1 > 0, log(freqs1 / freqs2), 0)
        KL = sum(freqs1 * LR)

        if (unit == "log2")
          KL = KL / log(2)  # change from log to log2 scale
        if (unit == "log10")
          KL = KL / log(10) # change from log to log10 scale

        return(KL)
      }

    JSdist <-
      function(x, y){
        0.5 * KLdist(x, (x + y) / 2) + 0.5 * KLdist(y, (x + y) / 2)
      }


    for (x in 2:(ncol(DATA.DF) - 2)) {
      DATA.DF[, x] <- as.numeric(as.character(DATA.DF[, x]))
    }

    # remove duplicate gene reads
    if (MULTI_PROBES == "highest") {
      temp1.df <- DATA.DF
      temp1.df$rowsum <- rowSums(temp1.df[, 2:(ncol(temp1.df) - 2)])
      temp1.df <-
        temp1.df[order(temp1.df$rowsum,
                       decreasing = TRUE,
                       na.last = NA),]
      row.names(temp1.df) <- 1:nrow(temp1.df)
      temp1.df <-
        (temp1.df[!duplicated(temp1.df[, (ncol(temp1.df) - 2)]),])
      temp1.df$rowsum <- NULL
      DATA.DF <-
        temp1.df[, c((ncol(temp1.df) - 1), 2:(ncol(temp1.df) - 2), (ncol(temp1.df)))]
    }
    else if (MULTI_PROBES == "mean") {
      DATA.DF[, 1] <- NULL
      keys <- colnames(DATA.DF)[(ncol(DATA.DF) - 1)]
      temp.dt <- data.table::as.data.table(DATA.DF)
      temp2 <- temp.dt[, lapply(data.table::.SD, mean), by = keys]
      DATA.DF <- as.data.frame(temp2)
    }
    else{
      message("Please choose to either keep the highest probe expression values or average them out")
      # terminate package
    }

    # format of data after above
    #  1st col #####  many cols ####    last col ###
    #  symbol    |     expression      | length

    if (FILTERED) {
      isexpr <-
        rowSums(edgeR::cpm(DATA.DF[2:(ncol(DATA.DF) - 1)]) > 1) >= 4
      DATA.DF <- DATA.DF[isexpr,]
    }

    # Quantile normalize if data is not NORMALIZED.
    if (NORMALIZED) {
      DATA.DF[2:(ncol(DATA.DF) - 1)] <-
        preprocessCore::normalize.quantiles(as.matrix(DATA.DF[2:(ncol(DATA.DF) - 1)]))
    }

    # sort by length
    DATA.DF$length <- as.numeric(as.character(DATA.DF$length))
    # sort dataframe by length in increasing order. also should remove the rows that had no length
    DATA.DF <-
      DATA.DF[order(DATA.DF$length, decreasing = FALSE, na.last = NA),]
    # redo the row names to match new order
    row.names(DATA.DF) <- 1:nrow(DATA.DF)
    # calculate number of steps needed for the data
    s <- ((nrow(DATA.DF) - WINDOW_SIZE) / STEP_SIZE) + 1
    #create dataframe for results
    simp.df <- data.frame("kb_length" = rep(0, s))
    # get labels of expression reads
    labels <-
      c("kb_length", colnames(DATA.DF)[2:(ncol(DATA.DF) - 1)])
    simp.df[labels] <- 0
    # create dataframe for pvalues
    pvalues.df <- data.frame("kb_length" = rep(0, s))
    pvalues.df[labels] <- 0

    # create dataframe for kL values
    JSdistances.df <- data.frame("kb_length" = rep(0, s))
    JSdistances.df[labels] <- 0


    # w is the counter
    w <- 1
    step_number <- 1
    while (step_number < s) {
      # get average length
      count <- 0
      if (WINDOW_STYLE == "mean") {
        count <- sum(DATA.DF$length[(w:(w + WINDOW_SIZE - 1))])
        count <- count / WINDOW_SIZE
      }
      if (WINDOW_STYLE == "median") {
        count <- median(DATA.DF$length[(w:(w + WINDOW_SIZE - 1))])
      }
      simp.df[step_number, 1] <- count
      JSdistances.df[step_number, 1] <- count
      pvalues.df[step_number, 1] <-  count
      for (column_index in (2:(ncol(DATA.DF) - 1))) {
        count <- 0
        if (WINDOW_STYLE == "mean") {
          count <- sum(DATA.DF[w:(w + WINDOW_SIZE - 1), column_index])
          count <- count / WINDOW_SIZE
        }
        if (WINDOW_STYLE == "median") {
          count <- median(DATA.DF[w:(w + WINDOW_SIZE - 1), column_index])
        }
        simp.df[step_number, column_index] <- count

        x <- DATA.DF[w:(w + WINDOW_SIZE - 1), CONTROL_COLUMN_INDEX]
        y <- DATA.DF[w:(w + WINDOW_SIZE - 1), column_index]
        pvalue <- stats::wilcox.test(x, y, exact = FALSE)$p.value
        if (!is.nan(pvalue)) {
          pvalues.df[step_number, column_index] <- pvalue
        } else {
          pvalues.df[step_number, column_index] <- 1
        }


      }
      w <- w + STEP_SIZE
      step_number <- step_number + 1
    }
    write.table(
      simp.df,
      "LONGO_out_final_table.tsv",
      sep = "\t",
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE
    )
    for(step_number in (2:s)){
    for (column_index in (2:(ncol(DATA.DF) - 1))) {
      x <- simp.df[2:step_number, CONTROL_COLUMN_INDEX]
      y <- simp.df[2:step_number, column_index]

      if (all(x == 0) && all(y == 0)) {
        JSdistance = 0
      } else {
        JSdistance <- JSdist(x, y)
      }

      if (is.na(JSdistance)) {
        JSdistance = 0
      }

      JSdistances.df[step_number, column_index] <- JSdistance
    }
    }
    # write.table(pvalues.df, "LONGO_out_pvalue_table.tsv", sep = "\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    return(list(simp.df, pvalues.df, JSdistances.df))
  }
