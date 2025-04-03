

# Clean the environment
rm(list = ls())

# Set the working directory
setwd("D:/Partners HealthCare Dropbox/Shi Sun/Matlab scripts/Pezaris/transcriptomy/human_lgn")

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(dplyr)
library(biomaRt)  # Add biomaRt for gene ID conversion
library(readxl)


# Load the processed Seurat object
data <- readRDS("data/lgn_human_brain_atlas.rds")


lgn_cells <- subset(data, subset = supercluster_term == "Thalamic excitatory")

# get gene names in the dataset
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- rownames(lgn_cells@assays$RNA@data)
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# only include the cosmx genes here
file_path <- "data/PMR-11548-Human-6k-Discovery-Gene-List.xlsx"
gene_data <- read_excel(file_path, sheet = "Gene and Probe Details", skip = 1)
gene_names <- gene_data$`Display Name`

# filter by genes present
keepgen <- gene_info$ensembl_gene_id[gene_info$hgnc_symbol %in% gene_names]

# filter out the dataset
lgn_cells <- subset(lgn_cells, features = keepgen)

# custom settings
ndims <- c(5, 10, 15, 20, 30) # pca ndims
feature <- c(50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000)
plot_list <- list()

count = 0
# loop between each set of cells
for (ndim in ndims){
  for (nfeatures in feature) {
    # lgn_cells <- NormalizeData(object = lgn_cells)
    lgn_cells <- FindVariableFeatures(object = lgn_cells, selection.method = "vst", nfeatures = nfeatures)
    lgn_cells <- ScaleData(object = lgn_cells, verbose = FALSE)
    lgn_cells <- RunPCA(object = lgn_cells, npcs = 30, verbose = FALSE)
    ElbowPlot(lgn_cells, ndims = 30)
    
    lgn_cells <- RunUMAP(object = lgn_cells, reduction = "pca", dims = 1:ndim)
    lgn_cells <- FindNeighbors(lgn_cells, reduction = "pca", dims = 1:ndim)
    lgn_cells <- FindClusters(lgn_cells, resolution = 0.2)
    
    umap_plot <- DimPlot(lgn_cells, reduction = "umap", group.by = "ident") + 
      ggtitle(paste0('pcas = ', ndim,', ngenes = ', nfeatures)) + 
      theme(aspect.ratio = 1,
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 7))
    
    count = count + 1
    plot_list[[count]] = umap_plot
    
  }
}


pdfname <- paste0("analysis_output/umap_iterations.pdf")
pdf(pdfname, height = 8, width = 15)
grid.arrange(grobs = plot_list, ncol = length(feature))

dev.off()

