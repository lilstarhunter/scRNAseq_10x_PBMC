# Data on local drive, no HPC 

library("Matrix")
library("dplyr")
library("Seurat")

# Load the PBMC dataset
# local laptop data directory = /Users/laurenstein/Desktop/local_data/scRNAseq_10x_PBMC/
pbmc.data <- Read10X(data.dir = "/Users/laurenstein/Desktop/local_data/scRNAseq_10x_PBMC/data/filtered_gene_bc_matrices/hg19")


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
metadata <- pbmc@meta.data


# Add number of genes per UMI for each cell to metadata
# Create new column called log10GenesPerUMI
pbmc$log10GenesPerUMI <- log10(pbmc$nFeature_RNA) / log10(pbmc$nCount_RNA)
pbmc$log10GenesPerUMI <- pbmc@meta.data$log10GenesPerUMI

# Add mitochondrial feature set
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Add mitoRatio to the meta.data table
pbmc$mitoRatio <- pbmc@meta.data$percent.mt / 100
pbmc@meta.data

# Create .RData object to load at any time
saveRDS(pbmc, file="/Users/laurenstein/Desktop/local_data/scRNAseq_10x_PBMC/data/pbmc_filtered_seurat.rds")



