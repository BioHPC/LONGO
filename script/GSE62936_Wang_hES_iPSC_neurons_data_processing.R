# Wang et al., 2015
# Raw data obtained from: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62936/matrix/

########################################################################################
# Raw data processing
########################################################################################

# Process the raw data so that all FPKM are in one file:

temp.1.df <- read.table(gzfile("GSM1536664_LW01_genes.csv.gz"), sep = ",", header = TRUE)
temp.2.df <- read.table(gzfile("GSM1536665_LW02_genes.csv.gz"), sep = ",", header = TRUE)
temp.3.df <- read.table(gzfile("GSM1536666_LW03_genes.csv.gz"), sep = ",", header = TRUE)
temp.4.df <- read.table(gzfile("GSM1536667_LW04_genes.csv.gz"), sep = ",", header = TRUE)
temp.5.df <- read.table(gzfile("GSM1536668_LW05_genes.csv.gz"), sep = ",", header = TRUE)
temp.6.df <- read.table(gzfile("GSM1536669_LW06_genes.csv.gz"), sep = ",", header = TRUE)
temp.7.df <- read.table(gzfile("GSM1536670_LW07_genes.csv.gz"), sep = ",", header = TRUE)
temp.9.df <- read.table(gzfile("GSM1536671_LW09_genes.csv.gz"), sep = ",", header = TRUE)
temp.10.df <- read.table(gzfile("GSM1536672_LW10_genes.csv.gz"), sep = ",", header = TRUE)
temp.11.df <- read.table(gzfile("GSM1536673_LW11_genes.csv.gz"), sep = ",", header = TRUE)
temp.12.df <- read.table(gzfile("GSM1536674_LW12_genes.csv.gz"), sep = ",", header = TRUE)
temp.13.df <- read.table(gzfile("GSM1536675_LW13_genes.csv.gz"), sep = ",", header = TRUE)
temp.14.df <- read.table(gzfile("GSM1536676_LW14_genes.csv.gz"), sep = ",", header = TRUE)

temp.1.df <- cbind(as.character(temp.1.df$gene_id), temp.1.df$FPKM)
temp.2.df <- cbind(as.character(temp.2.df$gene_id), temp.2.df$FPKM)
temp.3.df <- cbind(as.character(temp.3.df$gene_id), temp.3.df$FPKM)
temp.4.df <- cbind(as.character(temp.4.df$gene_id), temp.4.df$FPKM)
temp.5.df <- cbind(as.character(temp.5.df$gene_id), temp.5.df$FPKM)
temp.6.df <- cbind(as.character(temp.6.df$gene_id), temp.6.df$FPKM)
temp.7.df <- cbind(as.character(temp.7.df$gene_id), temp.7.df$FPKM)
temp.9.df <- cbind(as.character(temp.9.df$gene_id), temp.9.df$FPKM)
temp.10.df <- cbind(as.character(temp.10.df$gene_id), temp.10.df$FPKM)
temp.11.df <- cbind(as.character(temp.11.df$gene_id), temp.11.df$FPKM)
temp.12.df <- cbind(as.character(temp.12.df$gene_id), temp.12.df$FPKM)
temp.13.df <- cbind(as.character(temp.13.df$gene_id), temp.13.df$FPKM)
temp.14.df <- cbind(as.character(temp.14.df$gene_id), temp.14.df$FPKM)

data.df <- cbind(as.character(temp.1.df$gene_id), temp.1.df$FPKM,
                 temp.2.df$FPKM,
                 temp.3.df$FPKM,
                 temp.4.df$FPKM,
                 temp.5.df$FPKM,
                 temp.6.df$FPKM,
                 temp.7.df$FPKM,
                 temp.9.df$FPKM,
                 temp.10.df$FPKM,
                 temp.11.df$FPKM,
                 temp.12.df$FPKM,
                 temp.13.df$FPKM,
                 temp.14.df$FPKM
                 )
data.df <- as.data.frame(data.df)
colnames(data.df) <- c("gene_symbol",
                       "1",
                       "2",
                       "3",
                       "4",
                       "5",
                       "6",
                       "7",
                       "9",
                       "10",
                       "11",
                       "12",
                       "13",
                       "14")


write.table(data.df, file = "GSE62936_Wang_hES_iPSC_neurons_FPKM.csv", sep = ",") # Save file output

########################################################################################
# LONGO processing information
########################################################################################

# Species: Hsapiens
# Gene ID: external_gene_name
# Quantile normalized
# Data filtered (>1 cpm for at least 4 samples)
# Sliding median
# Multi probe values averaged
# Bin size: 200
# Step size: 40
