#' LONGO call biomaRt
#'
#' This is a internal function. This function calls the biomaRt database to
#' retrieve the gene symbols and lengths.
#'
#' @param data.df dataframe that contains the probe ids that will be checked
#' @param library_type the gene identifier of the data based off the
#' biomaRt database
#' @param ensembl BiomaRt mart based on the given species
#' @importFrom biomaRt getBM
#' @return returns the data.df with two new columns, length and symbol
callbiomaRt <- function(data.df, library_type, ensembl) {
  #call biomart for gene symbols and lengths from user identifier
  source_data.bm <-
    biomaRt::getBM(
      attributes = c(
        library_type,
        "external_gene_name",
        "start_position",
        "end_position"
      ),
      filters = library_type,
      values = data.df[, 1],
      mart = ensembl
    )
  #make sure the values are numbers
  source_data.bm[, 3] <-
    as.numeric(as.character(source_data.bm[, 3]))
  source_data.bm[, 4] <-
    as.numeric(as.character(source_data.bm[, 4]))

  #copy identifier, gene symbol
  affy.df <- source_data.bm[, 1:2]

  #add lengths
  #divide by 1000 to get kb
  affy.df$KBlength <-
    ((source_data.bm[, 4] - source_data.bm[, 3] + 1) / 1000)

  #order by length remove duplicate symbols
  #keeping the one that has the largest length
  affy.df <-
    affy.df[order(affy.df$KBlength,
                  decreasing = TRUE,
                  na.last = NA),]
  row.names(affy.df) <- 1:nrow(affy.df)
  affy.df <- (affy.df[!duplicated(affy.df[, 1]),])
  row.names(affy.df) <- 1:nrow(affy.df)

  return(affy.df)
}
