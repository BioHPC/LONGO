#' File opener
#'
#' This is an internal functuon. This function opens a file
#' and reads the data into a dataframe.
#'
#' @param fileLocation location to the file to be opened
#' @param separator variable that is used to separate columns
#' @param header true/false used to open a file with or without a header
#' @param commentChar character used as commented lines to be ignored
#' in the file
#' @return Returns nothing
openfile <- function(fileLocation, separator, header, commentChar) {
    data.df <-
        read.csv(
            file = fileLocation,
            header = header,
            sep = separator,
            blank.lines.skip = TRUE,
            comment.char = commentChar,
            na.strings = c("NA", " ", "")
        )
}
