#' openfile
#'
#'opens a file and reads the data into a dataframe
#' @param FILE location to the file to be opened
#' @param SEP variable that is used to separate columns
#' @param HEADER true/false used to open a file with or without a header
#' @param COMMENT_CHAR character used as commented lines to be ignored in the file
openfile <- function(FILE, SEP, HEADER, COMMENT_CHAR) {
  #add error checking to confirm it is opened
  #data.df <- read.table(FILE, sep=SEP, header=TRUE , blank.lines.skip=TRUE, comment.char="!", fill=TRUE)
  data.df <-
    read.csv(
      file = FILE,
      header = HEADER,
      sep = SEP,
      comment.char = COMMENT_CHAR,
      na.strings = c("NA", " ", "")
    )
}
