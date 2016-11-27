# LONGO: Gene Length-Dependent Expression Analysis Tool in Neuronal Cells

# Pre-requisites:

*   biomaRt package
    *   install by using the following R commands:  
            > source("https://bioconductor.org/biocLite.R")  
            > biocLite("biomaRt")
*   edgeR package
    *   install by using the following R commands:  
            > source("https://bioconductor.org/biocLite.R")  
            > biocLite("edgeR"")
*   preprocessCore package
    *   install by using the following R commands:  
            > source("https://bioconductor.org/biocLite.R")  
            > biocLite("preprocessCore")
*   shiny package:
    *   install by using the following R commands:  
            > install.packages("shiny")  
*   DT package
    *   install by using the following R commands:  
            > install.packages("DT")  
*   data.table package
    *   install by using the following R commands:  
            > install.packages("data.table")

*   hash package:
    *   install by using the following R command:  
            > install.packages("hash")

#Usage with shiny:
*   Launch the interface  
        > LONGO()
*   Set the initial variables
    *   Species gene ensembl: The species database used with the data file
    *   Gene identifier: Identifier in the data used to find gene
        names and lengths
    *   File: Data file with the first column being the gene identifier
        from above
    *   Header: True/False option
    *   Separator: Choose how the data file is structured
    *   Normalize: True/False option
    *   Filter RNAseq: True/False option
*   Submit data to be analyzed
    *   Process can take about a minute to communicate with the biomaRt 
        database
*   After initial analysis
    *   Alter analysis variables to get different results

#Usage without shiny:
*   Use the following R command:  
        > LONGOcmd(fileLocation = path_to_file, {separator = ","},
            {header = TRUE}, {commentChar = "!"},
            {species = "hsapiens_gene_ensembl"}, 
            {libraryType = "affy_primeview"}, {multiProbes = "highest"},
            {windowSize = 200}, {stepSize = 40}, {windowStyle = "mean"},
            {filterData = TRUE}, {normalizeData = TRUE}, {controlColumn = 2})
*   Output files are written to the working directory
