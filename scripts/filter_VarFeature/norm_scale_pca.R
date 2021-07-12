library(future)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ape)

options(warn = -1) #Warning about underscored features converted to dashes ignored

input = "/home/steinlm/scRNAseq_10x_PBMC/data/pbmc_ser_filtered_v1.rds"
output = "/home/steinlm/scRNAseq_10x_PBMC/data/pbmc_ser_filtered_tree_v1.rds"
image_output = "/home/steinlm/scRNAseq_10x_PBMC/plots/cluster/"
filter_df <- readRDS(input)

#Check and coerce seurat object into a dataframe
all_genes <- as.data.frame(rownames(filter_df))

#Return only non-mitochondrial genes
colnames(all_genes) <- c("genes")
nomtgenes <- all_genes[!grepl("^MT-",all_genes$genes),]
filter_df <- subset(filter_df, features = nomtgenes)

# head(filter_df)

plan("multiprocess",workers=16) #multiprocess is deprecated
options(future.globals.maxSize = 1000000 * 1024^2)

#Normalize and scale the data
filter_df <- NormalizeData(filter_df, normalization.method = "LogNormalize", scale.factor = 10000)
filter_df <- ScaleData(filter_df, features = nomtgenes, vars.to.regress = c("nCount_RNA"))
#
plan("sequential")

# Find Variable Features
filter_df <- FindVariableFeatures(filter_df, selection.method = "mvp", mean.cutoff=c(0.003,2))
var <- VariableFeatures(filter_df)
length(var)


png(filename = "/home/steinlm/scRNAseq_10x_PBMC/plots/cluster/filter_df_varFeat.png")
VariableFeaturePlot(object = filter_df)
dev.off()

# Clustering
filter_df <- RunPCA(filter_df, features = VariableFeatures(object = filter_df), ndims.print = 1:20, nFeatures.print = 20, npcs = 100)

png(filename = "/home/steinlm/scRNAseq_10x_PBMC/plots/cluster/filter_df_Elbow.png")
ElbowPlot(filter_df, ndims = 100)
dev.off()
#

saveRDS(filter_df, file = output)


filter_df <- FindNeighbors(filter_df, dims = 1:20)

filter_df <- FindClusters(filter_df, resolution = 0.5)

filter_df <- BuildClusterTree(filter_df, reorder = TRUE, reorder.numeric = TRUE, dims = 1:20)

Idents(filter_df) <- 'tree.ident'

levels(filter_df) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
# 
png("/home/steinlm/scRNAseq_10x_PBMC/plots/cluster/filter_df_clusterTree.png")
PlotClusterTree(filter_df)
dev.off()
# 
png("/home/steinlm/scRNAseq_10x_PBMC/plots/cluster/filter_df_VlnTree.png",  width = 720)
VlnPlot(filter_df, features = c("nCount_RNA"), pt.size=0, group.by='tree.ident') + theme(legend.position = 'none') #Grouped by cluster with no legend
dev.off()



# table(Idents(filter_df), filter_df@meta.data$subject) #Specify sample
# table(Idents(filter_df),filter_df@meta.data$group) #Specify treatment group

features = c("GFRAL", "GDF15", "FCN1", "LGALS2","PF4", "PPBP", "GZMB", "PRF1","GIMAP7", "AQP3", "GZMK", "LYAR","ABI3", "C1QA","PLBD1", "SELL", "TMEM222", "LBR","MAEA", "ATF7IP2","PHF5A", "MAPK1IP1L")
# 
p <- DotPlot(filter_df, features=features, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
  guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
  theme(
    axis.text.y=element_text(hjust=0,face="bold"),
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  )

ggsave(p, filename="/home/steinlm/scRNAseq_10x_PBMC/plots/cluster/dotplot_full1.png")

