library(future)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ape)

input = "/home/steinlm/scRNAseq_10x_PBMC/data/pbmc_ser_filtered_v1.rds"
output = "/home/steinlm/scRNAseq_10x_PBMC/data/pbmc_ser_filtered_tree_v1.rds"
image_output = "/home/steinlm/scRNAseq_10x_PBMC/plots/cluster"
filter_df <- readRDS(input)

#Check and coerce seurat object into a dataframe
all_genes <- as.data.frame(rownames(filter_df))

#Return only non-mitochondrial genes
colnames(all_genes) <- c("genes")
nomtgenes <- all.genes[!grepl("^MT-",all.genes$genes),]
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


png(paste(image_output,"filter_df_varFeat.png",sep=""))
VariableFeaturePlot(object = filter_df)
dev.off()

# Clustering
filter_df <- RunPCA(hayesfull, features = VariableFeatures(object = filter_df), ndims.print = 1:10, nFeatures.print = 10, npcs = 100)
filter_df
# 
png(paste(image_output,"filter_df_Elbow.png",sep=""))
ElbowPlot(filter_df2, ndims = 100)
dev.off()
# 

saveRDS(filter_df, file = output)
# 

# 
# filter_df <- FindNeighbors(filter_df, dims = 1:10)
# 
# filter_df <- FindClusters(filter_df, resolution = 0.5)
# 
# filter_df <- BuildClusterTree(filter_df, reorder = TRUE, reorder.numeric = TRUE, dims = 1:10)
# 
# Idents(filter_df) <- 'tree.ident'
# 
# levels(hayesfull) <- c(1,2,3,4,5,6,7,8,9,10)
# 
# png(paste(image_output,"filter_df_clusterTree.png",sep=""))
# PlotClusterTree(fitler_df)
# dev.off()
# 
# png(paste(image_output,"filter_df_VlnTree.png",sep=""),  width = 720)
# VlnPlot(filter_df, features = c("nCount_RNA"), pt.size=0, group.by='tree.ident') + theme(legend.position = 'none') #Grouped by cluster with no legend
# dev.off()



# table(Idents(filter_df), filter_df@meta.data$subject) #Specify sample
# table(Idents(filter_df),filter_df@meta.data$group) #Specify treatment group

# features = c()
# 
# p <- DotPlot(hayesfull, features=features, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
#   guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
#   theme(
#     axis.text.y=element_text(hjust=0,face="bold"),
#     axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
#     axis.title.x=element_blank(),
#     axis.title.y=element_blank()
#   )
# 
# ggsave(p, filename="/home/crist/Hayes_10x/dotplot_hayesfull1.png")
# 
# p <- DotPlot(hayesfull, features=features2, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
#   guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
#   theme(
#     axis.text.y=element_text(hjust=0,face="bold"),
#     axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
#     axis.title.x=element_blank(),
#     axis.title.y=element_blank()
#   )
# 
# ggsave(p, filename="/home/crist/Hayes_10x/dotplot_hayesfull2.png")
