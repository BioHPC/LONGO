# Tasic et al., 2016 Brain-map
# Raw data obtained from: http://casestudies.brain-map.org/celltax/download_archive

########################################################################################
# Raw data pre-processing
########################################################################################

data.df <- read.table(file = "genes_counts.csv", sep = ",", header = TRUE, row.names = 1)
meta.df <- read.table(file = "cell_metadata.csv", sep = ",", header = TRUE)
meta.df <- meta.df[,c(1,12)]
data.df <- as.data.frame(t(data.df))
data.df$long_name <- rownames(data.df)
merged.df <- merge(meta.df, data.df, "long_name")
data.df <- merged.df[,-1]
aggregated.df <- aggregate(data.df, by = list(data.df[,1]), FUN = function(x) {mean(as.numeric(as.character(x)))})
data.mat <- t(aggregated.df)
colnames(data.mat) <- data.mat[1,]
data.mat <- data.mat[-c(1:2),]
data.df <- as.data.frame(data.mat)
write.table(data.df, file = "Brain_atlas_mouse_visual_cortex_genes_RPKM_SubtypeMeans.csv", sep = ",")

########################################################################################
# Raw data processing
########################################################################################

# Process the raw data so that samples are averaged according to their condition

data.filename <- "Brain_atlas_mouse_visual_cortex_genes_RPKM_SubtypeMeans.csv"
out.filename <- unlist(strsplit(data.filename, "[.]"))[1] # Store filename to create output file later
data.df <- read.table(gzfile(data.filename), sep=",", header=TRUE , blank.lines.skip=TRUE)
data.df$Gene_symbol <- rownames(data.df)
data.df <- cbind(data.df[,ncol(data.df)], data.df[,-c(ncol(data.df))])
colnames(data.df)[1] <- "Gene_symbol"
write.table(data.df, file = paste(out.filename, "processed.csv", sep = "_"), sep = ",") # Save file output

########################################################################################
# LONGO processing information
########################################################################################

# Quantile normalized
# Data filtered (>1 cpm for at least 4 samples)
# Sliding median
# Multi probe values averaged
# Bin size: 200
# Step size: 40
# Control column: Astrocytes
