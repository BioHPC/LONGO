# Yu et al., Nature Communications 2014
# Data obtained from: http://pgx.fudan.edu.cn/ratbodymap/ratbodymap_data/fordownload/Ratbodymap_Gene_FPKM_v2.txt.gz

########################################################################################
# Raw data processing
########################################################################################

# Process the raw data so that samples are averaged according to their tissues, regardless of age.
# Note: the age does not affect the outcome of this particular analysis.

data.filename <- "Ratbodymap_Gene_FPKM_v2.txt.gz"
out.filename <- unlist(strsplit(data.filename, "[.]"))[1] # Store filename to create output file later
data.df <- read.table(gzfile("Ratbodymap_Gene_FPKM_v2.txt.gz"), header = TRUE, sep = "\t")
data.mat <- t(data.df)
colnames(data.mat) <- data.mat[1,]
data.mat <- data.mat[-1,]
data.df <- as.data.frame(data.mat)
data.df$tissue <- substring(rownames(data.df), 1, 3)
# data.df$tissue <- paste(substring(rownames(data.df), 1, 3), substring(rownames(data.df), 7, 9), sep = "") # Altenratively, split by age as well. Has no effect on outcome though.
data.df <- aggregate(data.df, by = list(data.df$tissue), FUN = function(x) {mean(as.numeric(as.character(x)))})
data.df$tissue <- NULL
data.mat <- t(data.df)
colnames(data.mat) <- data.mat[1,]
data.mat <- data.mat[-1,]
data.df <- as.data.frame(data.mat)
data.df$Gene_symbol <- rownames(data.df)
data.df <- cbind(data.df[,ncol(data.df)], data.df[,-ncol(data.df)])
colnames(data.df)[1] <- "Gene_symbol"
write.table(data.df, file = paste(out.filename, "TissueMeans.csv", sep = "_"), sep = ",")

########################################################################################
# LONGO processing information
########################################################################################

# Quantile normalized
# Sliding median
# Multi probe values averaged
# Bin size: 200
# Step size: 40
# Control column: Thymus
