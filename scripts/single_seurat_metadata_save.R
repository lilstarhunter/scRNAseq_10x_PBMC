#Loading a Single Dataset

library(Matrix)
library(dplyr)
library(Seurat)


# ==== Load the DATASET ====== #
df.data <- Read10X(data.dir = "input.path_filename")
```
# Initialize the Seurat object with the raw (non-normalized data).
df <- CreateSeuratObject(counts = df.data, project = "name", min.cells = 3, min.features = 200)


# === Create MetaData ==== #
metadata <- df@meta.data
# Add number of genes per UMI for each cell to metadata
df$log10GenesPerUMI <- log10(df$nFeature_RNA) / log10(df$nCount_RNA)
# Create new column called log10GenesPerUMI
df$log10GenesPerUMI <- df@meta.data$log10GenesPerUMI

# Add mitochondrial feature set - Make sure pattern is based on species ex. mouse is ^mt
df$percent.mt <- PercentageFeatureSet(df, pattern = "^MT-")

#Add mitoRatio to the meta.data table
df$mitoRatio <- df@meta.data$percent.mt / 100

# Create .rds object to load at any time
saveRDS(df, file="output.path_filename")



