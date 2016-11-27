#' LONGO call biomaRt
#'
#' This is a internal function. This function calls the biomaRt database to
#' retrieve the gene symbols and lengths.
#'
#' @param data.df dataframe that contains the probe ids that will be checked
#' @param libraryType the gene identifier of the data based off the
#' biomaRt database
#' @param ensembl BiomaRt mart based on the given species
#' @importFrom biomaRt getBM
#' @return returns the data.df with two new columns, length and symbol
callbiomaRt <- function(data.df, libraryType, ensembl) {
    biomaRt.df <-
        biomaRt::getBM(
            #libraryType.
            attributes=c(libraryType, "external_gene_name",
                "start_position", "end_position"
            ),
            filters=libraryType,
            values=data.df[, 1],
            mart=ensembl
        )

    if(is.na(biomaRt.df[2,2]) | nrow(biomaRt.df)<200){
        return(NULL)
    }
    #make sure the values are numbers
    biomaRt.df[, 3] < -as.numeric(as.character(biomaRt.df[, 3]))
    biomaRt.df[, 4] <- as.numeric(as.character(biomaRt.df[, 4]))

    #copy identifier, gene symbol
    completeBiomart.df <- biomaRt.df[, 1:2]

    #add lengths
    #divide by 1000 to get kb
    completeBiomart.df$KBlength <-((biomaRt.df[, 4] - biomaRt.df[, 3] + 1)
                        / 1000)

    # order by length remove duplicate symbols
    # keeping the one that has the largest length
    completeBiomart.df <- completeBiomart.df[order(completeBiomart.df$KBlength,
        decreasing=TRUE, na.last=NA),]
    row.names(completeBiomart.df) <- 1:nrow(completeBiomart.df)
    completeBiomart.df <-
        (completeBiomart.df[!duplicated(completeBiomart.df[, 1]),])
    row.names(completeBiomart.df) <- 1:nrow(completeBiomart.df)
    return(completeBiomart.df)
}
