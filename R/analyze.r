#' analyze data
#' Takes in a dataframe analyzes it. internal function
#'  @param data.df dataframe to be analyzed
#'  @param multi_probes Option to select what to do with multiple probes to a single gene, default=highest
#'  @param window_size Option to alter the size of the windows to be created from the dataframe, default=200
#'  @param step_size Option to alter the size of the steps used to create the windows, default=40
#'  @param window_style Option to alter the method used to create windows, default=mean
#'  @param filtered Option to filter data, default=TRUE
#'  @param normalized Option to normalize data, default=TRUE
#'  @param control_column_index column of the data used as a control for the statistical analyses
#'  @importFrom edgeR cpm
#'  @importFrom preprocessCore normalize.quantiles
#'  @importFrom data.table as.data.table
#'  @importFrom data.table as.data.frame
#'  @importFrom data.table .SD
#'  @importFrom data.table duplicated
#'  @return returns the simplified data
analyze <-
  function(data.df,
           multi_probes,
           window_size,
           step_size,
           window_style,
           filtered,
           normalized,
           control_column_index) {
   # library(edgeR)
   # library(preprocessCore)
   # library(data.table)

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
      function(x, y)
        0.5 * KLdist(x, (x + y) / 2) + 0.5 * KLdist(y, (x + y) / 2)


    for (x in 2:(ncol(data.df) - 2)) {
      data.df[, x] <- as.numeric(as.character(data.df[, x]))
    }

    # remove duplicate gene reads
    if (multi_probes == "highest") {
      temp1.df <- data.df
      temp1.df$rowsum <- rowSums(temp1.df[, 2:(ncol(temp1.df) - 2)])
      temp1.df <-
        temp1.df[order(temp1.df$rowsum,
                       decreasing = TRUE,
                       na.last = NA),]
      row.names(temp1.df) <- 1:nrow(temp1.df)
      temp1.df <-
        (temp1.df[!data.table::duplicated(temp1.df[, (ncol(temp1.df) - 2)]),])
      temp1.df$rowsum <- NULL
      data.df <-
        temp1.df[, c((ncol(temp1.df) - 1), 2:(ncol(temp1.df) - 2), (ncol(temp1.df)))]
    }
    else if (multi_probes == "mean") {
      data.df[, 1] <- NULL
      keys <- colnames(data.df)[(ncol(data.df) - 1)]
      temp.dt <- data.table::as.data.table(data.df)
      temp2 <- temp.dt[, lapply(data.table::.SD, mean), by = keys]
      data.df <- data.table::as.data.frame(temp2)
    }
    else{
      message("Please choose to either keep the highest probe expression values or average them out")
      # terminate package
    }

    # format of data after above
    #  1st col #####  many cols ####    last col ###
    #  symbol    |     expression      | length

    if (filtered) {
      isexpr <-
        rowSums(edgeR::cpm(data.df[2:(ncol(data.df) - 1)]) > 1) >= 4
      data.df <- data.df[isexpr,]
    }

    # Quantile normalize if data is not normalized.
    if (normalized) {
      data.df[2:(ncol(data.df) - 1)] <-
        preprocessCore::normalize.quantiles(as.matrix(data.df[2:(ncol(data.df) - 1)]))
    }

    # sort by length
    data.df$length <- as.numeric(as.character(data.df$length))
    # sort dataframe by length in increasing order. also should remove the rows that had no length
    data.df <-
      data.df[order(data.df$length, decreasing = FALSE, na.last = NA),]
    # redo the row names to match new order
    row.names(data.df) <- 1:nrow(data.df)
    # calculate number of steps needed for the data
    s <- ((nrow(data.df) - window_size) / step_size) + 1
    #create dataframe for results
    simp.df <- data.frame("kb_length" = rep(0, s))
    # get labels of expression reads
    labels <-
      c("kb_length", colnames(data.df)[2:(ncol(data.df) - 1)])
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
      if (window_style == "mean") {
        count <- sum(data.df$length[(w:(w + window_size - 1))])
        count <- count / window_size
      }
      if (window_style == "median") {
        count <- median(data.df$length[(w:(w + window_size - 1))])
      }
      simp.df[step_number, 1] <- count
      JSdistances.df[step_number, 1] <- count
      pvalues.df[step_number, 1] <-  count
      for (column_index in (2:(ncol(data.df) - 1))) {
        count <- 0
        if (window_style == "mean") {
          count <- sum(data.df[w:(w + window_size - 1), column_index])
          count <- count / window_size
        }
        if (window_style == "median") {
          count <- median(data.df[w:(w + window_size - 1), column_index])
        }
        simp.df[step_number, column_index] <- count

        x <- data.df[w:(w + window_size - 1), control_column_index]
        y <- data.df[w:(w + window_size - 1), column_index]
        pvalue <- stats::wilcox.test(x, y, exact = FALSE)$p.value
        if (!is.nan(pvalue)) {
          pvalues.df[step_number, column_index] <- pvalue
        } else {
          pvalues.df[step_number, column_index] <- 1
        }

        if (all(x == 0) && all(y == 0)) {
          JSdistance = 0
        } else {
          JSdistance <- JSdist(x, y)
        }

        if (is.na(JSdistance)) {
          JSdistance = 0
        }

        JSdistances.df[step_number, column_index] <- JSdistance
        # pvalues.df[step_number,column_index] <- pvalue
      }
      w <- w + step_size
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

    # write.table(pvalues.df, "LONGO_out_pvalue_table.tsv", sep = "\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    return(list(simp.df, pvalues.df, JSdistances.df))
  }
