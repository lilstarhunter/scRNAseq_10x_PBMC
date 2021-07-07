# SCRIPT FUNCITON: read a 10x dataset, create a seurat object, and add metadata 
# OUTPUT: RDS file containing new metadata and a seurat object


# Dataset Specifications
# 1. Single dataset, single batch (no merging required)
# 2. TAR of RAW data with barcodes.tsv, matrix.mtx and features.tsv
# 3. No biological/experimental groups added at this point
# 4. Human Query
# *** NOTE: RDS files at this step will all be named with the suffix seraut_meta_v1.rds *** #


library(Seurat)
library(dplyr)
library(Matrix)


options(warn = -1) #Warning about underscored features converted to dashes ignored

# ============================ #
# ==== Load the DATASET ====== #
# ============================ #
# *** NOTE:  ALL FILES / PATHS ON PMACS HPC
path =  "/home/steinlm/scRNAseq_10x_PBMC/data/filtered_gene_bc_matrices/hg19"
df.data <- Read10X(data.dir = path)

# Initialize the Seurat object with the raw (non-normalized data).
df <- CreateSeuratObject(counts = df.data, project = "pbmc", min.cells = 3, min.features = 200) #Make sure to change name of project to something relevant


# ============================ #
# === Create MetaData ======== #
# ============================ #
metadata <- df@meta.data

# Add number of genes per UMI for each cell to metadata
df$log10GenesPerUMI <- log10(df$nFeature_RNA) / log10(df$nCount_RNA)

# Create new column called log10GenesPerUMI
df$log10GenesPerUMI <- df@meta.data$log10GenesPerUMI

# Add mitochondrial feature set - Make sure pattern is based on species ex. mouse is ^mt
df$percent.mt <- PercentageFeatureSet(df, pattern = "^MT-")

#Add mitoRatio to the meta.data table
df$mitoRatio <- df@meta.data$percent.mt / 100


# ================================ #
# ===== Save RDS/MetaData ======== #
# ================================ #
output_path = "/home/steinlm/scRNAseq_10x_PBMC/data/pbmc_ser_meta_v1.rds"
saveRDS(df, file=output_path)



