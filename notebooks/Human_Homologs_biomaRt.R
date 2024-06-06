# Install the biomaRt package if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

# Load the biomaRt library
library(biomaRt)

# Use Ensembl's bioMart service
ensembl <- useEnsembl(biomart = "ensembl", mirror = "useast")

# Select the dataset for yeast
yeast_dataset <- useDataset("scerevisiae_gene_ensembl", mart = ensembl)

# Select the dataset for human
human_dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Yeast genes list
yeast <- read.csv("file_path", header = FALSE)

#Pre-processing
yeast_1 <- yeast[2:nrow(yeast), 1:3]
colnames(yeast_1) <- yeast_1[1, ]
yeast_1 <- yeast_1[-1, ]

# Getting the gene list
yeast_genes<- yeast_1['Gene/ORF']

# Get orthologs
orthologs <- getLDS(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "external_gene_name",
                    values = yeast_genes,
                    mart = yeast_dataset,
                    attributesL = c("ensembl_gene_id", "external_gene_name"),
                    martL = human_dataset)

# View results
print(orthologs)
