#' dict
#'
#' creates and returns a dicitonary
#' @param data.df inputed data
#' @param source_data.bm biomaRt database
#' @importFrom hash keys
#' @importFrom hash .set
#' @return returns a dictionary
# @export
dict <- function(data.df, source_data.bm) {
  symbol.dict <- hash::hash(keys = data.df[, 1], values = NA)
  .set(symbol.dict, keys = source_data.bm[, 1], values = source_data.bm[, 2])
  len.dict <- hash::hash(keys = data.df[, 1], values = NA)
  .set(len.dict, keys = source_data.bm[, 1], values = source_data.bm[, 3])

  symdict.fun <- function(key) {
    return(symbol.dict[[key]])
  }
  lendict.fun <- function(key) {
    return(len.dict[[key]])
  }
  #get lengths and names
  data.df$symbol <- apply(as.matrix(data.df[, 1]), 1, symdict.fun)
  data.df$length <- apply(as.matrix(data.df[, 1]), 1, lendict.fun)
  #remove na's
  data.df <-
    data.df[!(is.na(data.df[, (ncol(data.df) - 1)]) |
                (data.df[, (ncol(data.df) - 1)] == "")), ]
  data.df <-
    data.df[!(is.na(data.df[, (ncol(data.df))]) |
                (data.df[, (ncol(data.df))] == "")), ]

  write.table(
    data.df,
    "LONGO_out_raw_table.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
  )
  return(data.df)
}

