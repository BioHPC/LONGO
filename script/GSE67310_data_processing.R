# Treutlein et al., Nature 2016
# Raw data obtained from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67310
# File name: GSE67310_iN_data_log2FPKM_annotated.txt.gz

########################################################################################
# Raw data processing
########################################################################################

# Process the raw data so that samples are averaged according to their condition

data.filename <- "GSE67310_iN_data_log2FPKM_annotated.txt.gz"
out.filename <- unlist(strsplit(data.filename, "[.]"))[1] # Store filename to create output file later
data.df <- read.table(gzfile(data.filename), sep="\t", header=TRUE , blank.lines.skip=TRUE)
data.df <- data.df[,-c(1,3:5)] # Remove extra information
data.df <- aggregate.data.frame(data.df, by = list(data.df$assignment), FUN = mean) # Average rows by sample condition
data.df <- data.df[,-2] # Remove extra information
data.mat <- t(data.df) # Switch columns to rows
colnames(data.mat) <- data.mat[1,] # Change column names to row 1
data.mat <- data.mat[-1,] # Remove extra information
data.df <- as.data.frame(data.mat) # Change back to dataframe
data.df <- cbind(rownames(data.df), data.df) # Make sure first column is gene symbols
colnames(data.df)[1] <- "Gene_symbol"
write.table(data.df, file = paste(out.filename, "SampleMeans.csv", sep = "_"), sep = ",") # Save file output

########################################################################################
# LONGO processing information
########################################################################################

# Quantile normalized
# Sliding median
# Multi probe values averaged
# Bin size: 200
# Step size: 40
# Control column: Fibroblast
