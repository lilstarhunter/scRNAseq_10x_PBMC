# Single-cell RNASeq Analysis - Peripheral Blood Mononuclear Cells

## 10X Genomics - Seurat

### Objective: Utilize the Seurat R Toolkit to analyze PBMC to identify distinct cell populations and explore differential expression

<img src="https://d2ygg2jwuhi4sz.cloudfront.net/wp/wp-content/uploads/2018/12/AdobeStock_208548494-1200x480.jpeg">

**[Data Source](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)**

## Step 1: Read Input Data

- `Read10x` function reads in the output of the 10x [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline

- Dependencies required: dplyr, Seurat, patchwork
