#' LONGO function
#'
#' This function allows the user to preset the variables and allows a faster
#' analysis of data. Please use the shiny interface first to get an
#' understanding of how this package analyzes data.
#'
#' @param FILE Address to the file to be analyzed
#' @param SEP Option to specify the separated value, default=","
#' @param HEADER Option to specify if there is a header in the file, default=TRUE
#' @param COMMENT_CHAR Option to specify the commented character in the loaded file, default="!"
#' @param SPECIES species name of the data based off the biomaRt database
#' @param LIBRARY_TYPE the gene identifier of the data based off the biomaRt database
#' @param MULTI_PROBES Option to change the method for handling multiple probes going to the same gene, default=highest
#' @param WINDOW_SIZE Option to alter the size of the windows to be created from the dataframe, default=200
#' @param STEP_SIZE Option to alter the size of the steps used to create the windows, default=40
#' @param WINDOW_STYLE Option to alter the method used to create windows, default=mean
#' @param FILTERED Option to filter data, default=TRUE
#' @param NORMALIZED Option to normalize data, default=TRUE
#' @param CONTROL_COLUMN_INDEX column of the data used as a control for the statistical analysis
#' @examples
#' LONGOcmd(datafile)
#' LONGOcmd(datafile, WINDOW_SIZE=100, STEP_SIZE=20)
#' LONGOcmd(datafile, MULTI_PROBES="mean")
#' @export
LONGOcmd <- function(FILE,
                     SEP = ",",
                     HEADER = TRUE,
                     COMMENT_CHAR = "!",
                     SPECIES = "hsapiens_gene_ensembl",
                     LIBRARY_TYPE = "affy_primeview",
                     MULTI_PROBES = "highest",
                     WINDOW_SIZE = 200,
                     STEP_SIZE = 40,
                     WINDOW_STYLE = "mean",
                     FILTERED = TRUE,
                     NORMALIZED = TRUE,
                     CONTROL_COLUMN_INDEX = 2) {
  message("* Reading file")
  if (!is.data.frame(FILE)) {
    data.df <- openfile(FILE, SEP, HEADER, COMMENT_CHAR)
  }
  else{
    #dataframe Probe ID first column expression values all others
    data.df <- FILE
  }


  message("* Getting gene symbol and length from biomaRt")
  species_ensembl <-
    biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                     host = "www.ensembl.org",
                     dataset = SPECIES)
  biomaRt.df <- callbiomaRt(data.df, LIBRARY_TYPE, species_ensembl)
  #biomaRt.df <- read.table("LONGO_out_affy_table.tsv", sep="\t", header=TRUE)  # For fast testing

  # Use dict symbol and length
  message("* Hashing symbol and length")
  data.df <- dict(data.df, biomaRt.df)

  # Analyze it
  message("* Analyzing the data")
  temp.df <-
    analyze(
      data.df,
      MULTI_PROBES,
      WINDOW_SIZE,
      STEP_SIZE,
      WINDOW_STYLE,
      FILTERED,
      NORMALIZED,
      CONTROL_COLUMN_INDEX
    )
  simp.df <- as.data.frame(temp.df[1])
  write.table(
    simp.df,
    "LONGO_out_final_table.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
  )
  p <- as.data.frame(temp.df[2])
  write.table(
    p,
    "LONGO_out_p_values.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
  )
  js <- as.data.frame(temp.df[3])
  write.table(
    js,
    "LONGO_out_js_values.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
  )
  lq <- as.data.frame(temp.df[4])
  write.table(
    lq,
    "LONGO_out_LONGO_Quotient.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
  )
  # Plot graph
  plotLONGO(simp.df)
  message("* Finished LONGO")
}
