# SCRIPT FUNCITON: Filtering dataset based on initial QC conditions
# OUTPUT: new Seurat object


# Dataset Specifications
# Run script after initialQC.R
# Input is the updated _ser_meta_v1.rds

library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)


# Set paths
input_path = "/home/steinlm/scRNAseq_10x_PBMC/data/pbmc_ser_meta_v1.rds"
output_path = "/home/steinlm/scRNAseq_10x_PBMC/plots/QC/filtered/"
data_output_path = "/home/steinlm/scRNAseq_10x_PBMC/data/pbmc_ser_filtered_v1.rds"

df <- readRDS(input_path)

#Load metdata
metadata <- df@meta.data

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA,
                mitoRatio = mitoRatio)


# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_df <- subset(x = df, 
                        subset= (nCount_RNA >= 500) & 
                          (nFeature_RNA >= 250) & 
                          (nFeature_RNA < 2500) & 
                          (log10GenesPerUMI > 0.80) & 
                          (mitoRatio < 0.20))


# =============================== #
# ===== GENE LEVEL FILTER ======= #
# =============================== #

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_df, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_df <- CreateSeuratObject(filtered_counts, meta.data = filtered_df@meta.data)

# =============================== #
# ===== REASSES QC METRIC ======= #
# =============================== #

# Visualize the number UMIs/transcripts per cell
#Indicates how deeply sequenced the samples were
metadata %>% 
  ggplot(aes(x=nUMI)) + 
  geom_histogram(aes(y=..density..), color = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "green") + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

ggsave(paste(output_path,"UMIperCell.png",sep=""))

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(x=nGene)) + 
  geom_histogram(aes(y=..density..), color = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "green") + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

ggsave(paste(output_path,"GenesperCell.png",sep=""))

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(y=log10(nGene))) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

ggsave(paste(output_path,"GenesperCell_Boxplot.png",sep=""))

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(x=percent.mt)) + 
  geom_histogram(aes(y=..density..), color = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "green") + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

ggsave(paste(output_path,"MitoperCell.png",sep=""))

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI 
# Look into novelty score, reference states it should be about 0.8
metadata %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_histogram(aes(y=..density..), color = "black", fill = "white") +
  geom_density(alpha = 0.2, fill ="green") +
  theme_classic() +
  geom_vline(xintercept = 0.8)

ggsave(paste(output_path,"GenesperUMI.png",sep=""))

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "purple", high = "yellow") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

ggsave(paste(output_path,"GenesCorrUMI.png",sep=""))

# =============================== #
# ===== NEW FILTERED RDS ======== #
# =============================== #
saveRDS(filtered_df, file=data_output_path)
