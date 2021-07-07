# SCRIPT FUNCITON: Create QC plots to determine initial filtering 
# OUTPUT: .PNG files from unfiltered dataset


# Dataset Specifications
# Run script after initial single_seurat_metadata_save.R
# Input is the generated RDS file


library(Matrix)
library(dplyr)
library(Seurat)
library(ggplot2)

# Set paths
input_path = "/home/steinlm/scRNAseq_10x_PBMC/data/pbmc_ser_meta_v1.rds"
output_path = "/home/steinlm/scRNAseq_10x_PBMC/plots/QC/unfiltered/"
df <- readRDS(input_path)

#Load metdata
metadata <- df@meta.data


# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA,
                mitoRatio = mitoRatio)

## =========================================== ##
## =========== UNIVARIATE QC PLOTS =========== ##
## =========================================== ##

# Visualize the number UMIs/transcripts per cell = Sequencing depth
metadata %>% 
  ggplot(aes(x=nUMI)) + 
  geom_histogram(aes(y=..density..), color = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "purple") + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

ggsave(paste(output_path,"UMIperCell.png",sep=""))

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(x=nGene)) + 
  geom_histogram(aes(y=..density..), color = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "purple") + 
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
  geom_density(alpha = 0.2, fill = "purple") + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

ggsave(paste(output_path,"MitoperCell.png",sep=""))

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI 
# Look into novelty score, reference states it should be about 0.8
metadata %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

ggsave(paste(output_path,"GenesperUMI.png",sep=""))

## =========================================== ##
## =========== MULTIVARIATE QC PLOTS =========== ##
## =========================================== ##

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

ggsave(paste(output_path,"GenesperUMI.png",sep=""))

