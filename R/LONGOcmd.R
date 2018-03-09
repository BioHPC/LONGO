#' LONGO function
#'
#' This function allows the user to preset the variables and allows a faster
#' analysis of data. Please use the shiny interface first to get an
#' understanding of how this package analyzes data.
#'
#' @param fileLocation Address to the file to be analyzed
#' @param separator Option to specify the separated value, default=","
#' @param header Option to specify if there is a header in the file,
#' default=TRUE
#' @param commentChar Option to specify the commented character in the loaded
#' file, default="!"
#' @param species species name of the data based off the biomaRt database
#' @param libraryType the gene identifier of the data based off the
#' biomaRt database
#' @param multiProbes Option to change the method for handling multiple
#' probes going to the same gene, default=mean
#' @param windowSize Option to alter the size of the windows to be created
#' from the dataframe, default=200
#' @param stepSize Option to alter the size of the steps used to create the
#' windows, default=40
#' @param windowStyle Option to alter the method used to create windows,
#' default=mean
#' @param filterData Option to filter data, default=TRUE
#' @param normalizeData Option to normalize data, default=TRUE
#' @param controlColumn column of the data used as a control for the
#' statistical analysis
#' @examples
#' LONGOcmd(exampleRatData, species = "mmusculus_gene_ensembl",
#' libraryType = "affy_mouse430a_2")
#' \donttest{
#' LONGOcmd("~/data.csv")
#' LONGOcmd("~/data.csv", windowSize=100, stepSize=20)
#' LONGOcmd("~/data.csv", multiProbes="highest")
#' }
#' @return Returns nothing
#' @export
LONGOcmd <- function(fileLocation,
                    separator=",",
                    header=TRUE,
                    commentChar="!",
                    species="hsapiens_gene_ensembl",
                    libraryType="affy_primeview",
                    multiProbes="mean",
                    windowSize=200,
                    stepSize=40,
                    windowStyle="mean",
                    filterData=TRUE,
                    normalizeData=TRUE,
                    controlColumn=2) {
    message("* Reading file")
    if (!is.data.frame(fileLocation)) {
        data.df <- openfile(fileLocation, separator, header, commentChar)
    }
    else{
        #dataframe Probe ID first column expression values all other cols
        data.df <- fileLocation
    }

    message("* Getting gene symbol and length from biomaRt")
    speciesEnsembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                        host="www.ensembl.org",
                                        dataset=species)
    biomaRt.df <- callbiomaRt(data.df, libraryType, speciesEnsembl)
    if(is.null(biomaRt.df)){
        message("The gene identifier or the species
            selected was incorrect for the uploaded data. Please select
            the correct species and gene identifier for your data.")
        return()
    }
    # Use dict symbol and length
    message("* Hashing symbol and length")
    data.df <- dict(data.df, biomaRt.df)

    # Analyze it
    message("* Analyzing the data")
    temp.df <- analyze(data.df, multiProbes, windowSize, stepSize,
                        windowStyle, filterData, normalizeData,
                        controlColumn)
    if(is.null(temp.df)){
        message("There were less than 200 genes identified
                with the species and gene identifier. Please make sure you have
                the correct inputs for your data and resubmit.")
        return()
    }

    simp.df <- as.data.frame(temp.df[1])
    write.table(simp.df, "LONGO_out_final_table.tsv", sep="\t",
                quote=FALSE, col.names=TRUE, row.names=FALSE)
    p <- as.data.frame(temp.df[2])
    write.table(p, "LONGO_Out_P_Values.tsv", sep="\t", quote=FALSE,
                col.names=TRUE,row.names=FALSE)
    js <- as.data.frame(temp.df[3])
    write.table(js, "LONGO_Out_Divergent_Values.tsv", sep="\t", quote=FALSE,
                col.names=TRUE, row.names=FALSE)
    lq <- as.data.frame(temp.df[4])
    write.table(lq, "LONGO_Out_Long_Gene_Quotient.tsv", sep="\t", quote=FALSE,
                col.names=TRUE, row.names=FALSE)
    # Plot graph
    plotLONGO(simp.df)
    message("* Finished LONGO")
}
