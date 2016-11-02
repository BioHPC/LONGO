# LONGO: Gene Length-Dependent Expression Analysis Tool in Neuronal Cells

# Pre-requisites
*    hash package:
    *   install by using the following R command:  
        \> install.packages("hash")
*    biomaRt package
    *   install by using the following R commands:  
        > source("https://bioconductor.org/biocLite.R")  
        > biocLite("biomaRt")


#Usage:
   LONGO(path_to_file, species, {MULTI_PROBES="highest"}, {WINDOW_SIZE=200}, {STEP_SIZE=40})

#Data file **must** have the following format:
* any notation about the probe reads must start with a "!"
* blank rows are ignored
* first non blank row must have the column names which will be used in the plot
* the first column must have the probe ID's
* all the following columns are reads

failure to follow the above requirments will cause the LONGO package to fail