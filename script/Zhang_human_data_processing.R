# Zhang et al., Neuron 2015
# Raw data obtained from: http://web.stanford.edu/group/barres_lab/brainseq2/brainseq2.html
# File name: TableS4-HumanMouseMasterFPKMList.xlsx

########################################################################################
# Raw data pre-processing
########################################################################################

# Gender and age information were removed from the .xlsx file, and the file was converted
# to .csv before editing here.

########################################################################################
# Raw data processing
########################################################################################

# Process the raw data so that samples are averaged according to their condition

#data.filename <- "GSE73721_Human_and_mouse_table.csv.gz"
data.filename <- "Zhang_Human_preprocessed.csv"
out.filename <- unlist(strsplit(data.filename, "[_p]"))[1] # Store filename to create output file later
data.df <- read.table(gzfile(data.filename), sep=",", header=FALSE , blank.lines.skip=TRUE)
data.mat <- t(data.df)
colnames(data.mat) <- data.mat[1,]
data.mat <- data.mat[-1,]
data.df <- as.data.frame(data.mat)
data.df <- aggregate.data.frame(data.df, by = list(data.df$Gene_symbol), FUN = function(x) {mean(as.numeric(as.character(x)))}) # Average rows by sample condition
data.df <- data.df[,-2] # Remove extra information
data.mat <- t(data.df) # Switch columns to rows
colnames(data.mat) <- data.mat[1,] # Change column names to row 1
data.mat <- data.mat[-1,] # Remove extra information
data.df <- as.data.frame(data.mat) # Change back to dataframe
data.df <- cbind(rownames(data.df), data.df) # Make sure first column is gene symbols
colnames(data.df)[1] <- "Gene_symbol"
data.df <- data.df[,-6] # Remove microglia/macrophage sample--the differences between microglia and other glial cell types is distracting
write.table(data.df, file = paste(out.filename, "SampleMeans.csv", sep = "_"), sep = ",") # Save file output


########################################################################################
# LONGO processing information
########################################################################################

# Species: Hspaiens
# Gene ID: external_gene_name
# Quantile normalized
# Data filtered (>1 cpm for at least 4 samples)
# Sliding median
# Multi probe values averaged
# Bin size: 200
# Step size: 40
