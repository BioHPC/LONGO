#' Calculates LONGO quotient
#'
#' This is an internal functuon. This function calculcates the LONGO quotient.
#'
#' @param data.df The simplified data from analysis
#' @param JSdata.df The JS distance data from analysis
#' @param controlColumn The column in the data which should be used as
#' control
#' @return Returns final data
LONGOquotient <- function(
    data.df,
    JSdata.df,
    controlColumn){

    maxJS <- max(JSdata.df[nrow(JSdata.df),2:ncol(JSdata.df)])
    for (i in 2:(ncol(JSdata.df))) {
        JSdata.df[,i] <- JSdata.df[,i]/maxJS
    }
    data.df <- data.df[,c(ncol(data.df),2:ncol(data.df))]

    JSdata.df[,controlColumn] <- 0
    LQdata.df <- JSdata.df[!(JSdata.df[,1]<150),]
    LQSignData <- data.df[!(JSdata.df[,1]<150),]
    LQSign <- data.frame(matrix(0, nrow = 1, ncol = ncol(data.df)))

    colnames(LQSign) <- colnames(LQdata.df)
    avgControl <- mean(LQSignData[,controlColumn])
    for (i in 2:ncol(LQSignData)) {
        avg <- mean(LQSignData[,i])
        if (avg >= avgControl) {
            LQSign[1,i] <- 1
        } else {
            LQSign[1,i] <- -1
        }
    }
    LQSign <- LQSign[,-1]
    for (i in 2:ncol(LQdata.df)) {
        min_LQ <- min(LQdata.df[,i])
        for (j in 1:nrow(LQdata.df)) {
            LQdata.df[j,i] <- LQdata.df[j,i] - min_LQ
        }
    }

    final.df <- data.frame(matrix(0, nrow = 1, ncol = ncol(LQdata.df)))
    for (i in 2:ncol(LQdata.df)) {
        final.df[1,i] <- LQdata.df[nrow(LQdata.df),i]-LQdata.df[1,i]
    }
    colnames(final.df) <- colnames(data.df)
    final.df <- final.df[,-1]
    for (i in 1:ncol(final.df)) {
        final.df[1,i] <- final.df[1,i]*LQSign[1,i]
    }
    final.df <- final.df[, order(as.numeric(final.df[1, ]))]
    return(final.df)
}
