# Miller et al., 2014 BrainSpan
# Subject: H376.IV.03, female, 21pcw
# Data source: http://www.brainspan.org/api/v2/well_known_file_download/278444094

########################################################################################
# Raw data processing
########################################################################################

# Process the raw data so that samples are averaged according to their structures
col_metadata <- read.table("columns_metadata.csv", header = TRUE, sep = ",")
col_metadata$structure_acronym <- gsub('[a-z]', '', col_metadata$structure_acronym) # remove lowercase letters
structures <- c("VZ", "SZ", "IZ", "SP", "CP", "MZ", "CGE", "LGE", "MGE", "SG") # Here are the structures we are interested in
col_metadata$col_number <- rownames(col_metadata) # Preserve row number so we can extract the relevant expression data
for (i in 1:nrow(col_metadata)) {
  if (col_metadata[i,3] %in% structures) {
    col_metadata[i,] <- col_metadata[i,]
  } else {
    col_metadata[i,] <- NA
  }
}
col_metadata <- na.omit(col_metadata)
data.df <- read.table("expression_matrix.csv", header = FALSE, sep = ",")
data.mat <- t(data.df)
data.df <- as.data.frame(data.mat)
data.df <- cbind(rownames(data.df), data.df)
colnames(data.df)[1] <- "col_number"
data.df$col_number <- gsub('V', '', data.df$col_number)
data.df[,1] <- as.numeric(data.df[,1]) -1
colnames(data.df) <- data.df[1,]
data.df <- data.df[-1,]
colnames(data.df)[1] <- "col_number"
data.df <- merge(col_metadata, data.df, by = "col_number")
data.df <- data.df[,-c(1:3,5)]
data.df <- aggregate.data.frame(data.df, by = list(data.df$structure_acronym), FUN = function(x) {mean(as.numeric(as.character(x)))})
data.df <- data.df[,-2]
data.mat <- t(data.df)
colnames(data.mat) <- data.mat[1,]
data.df <- as.data.frame(data.mat)
data.df <- data.df[-1,]
# Convert probe IDs to gene symbols
rows_metadata <- read.table("rows_metadata.csv", header = TRUE, sep = ",")
rows_metadata <- rows_metadata[,c(1,4)]
data.df$probeset_id <- rownames(data.df)
data.df <- merge(rows_metadata, data.df, by = "probeset_id")
data.df <- data.df[,-1]
# Write our processed data
write.table(data.df, file = "Brainspan_developing_human_microarray.processed.csv", sep = ",")


########################################################################################
# LONGO processing information
########################################################################################

# Quantile normalized
# Sliding median
# Multi probe values averaged
# Bin size: 200
# Step size: 40
# Control column: VZ
