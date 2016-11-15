#' File opener
#'
#' This is an internal functuon. This function opens a file
#' and reads the data into a dataframe.
#'
#' @param FILE location to the file to be opened
#' @param SEP variable that is used to separate columns
#' @param HEADER true/false used to open a file with or without a header
#' @param COMMENT_CHAR character used as commented lines to be ignored in the file
#' @return Returns nothing
openfile <- function(FILE, SEP, HEADER, COMMENT_CHAR) {
  data.df <-
    read.csv(
      file = FILE,
      header = HEADER,
      sep = SEP,
      comment.char = COMMENT_CHAR,
      na.strings = c("NA", " ", "")
    )
}
