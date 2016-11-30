# Mazzoni et al., 2011 Mouse ESC to motorneuron
# Source data: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19372/matrix/

########################################################################################
# Raw data processing
########################################################################################

data.df <- read.table(gzfile("GSE19372_series_matrix.txt.gz"), sep = "\t", header = TRUE, fill = TRUE, blank.lines.skip=TRUE, comment.char = "!")
write.table(data.df, file = "GSE19372_ESC2Motomin.csv", sep = ",")

########################################################################################
# LONGO processing information
########################################################################################

# Species: Mmusculus
# Gene ID: affy_mouse430_2
# Quantile normalized
# Sliding median
# Multi probe values averaged
# Bin size: 200
# Step size: 40
