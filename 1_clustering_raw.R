

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

# Load the processed Seurat object
data <- readRDS("data/lgn_human_brain_atlas.rds")


lgn_cells <- subset(data, subset = supercluster_term == "Thalamic excitatory")


# lgn_cells <- NormalizeData(object = lgn_cells)
lgn_cells <- FindVariableFeatures(object = lgn_cells, selection.method = "vst", nfeatures = 3000)
lgn_cells <- ScaleData(object = lgn_cells, verbose = FALSE)
lgn_cells <- RunPCA(object = lgn_cells, npcs = 30, verbose = FALSE)
ElbowPlot(lgn_cells)

lgn_cells <- RunUMAP(object = lgn_cells, reduction = "pca", dims = 1:15)
lgn_cells <- FindNeighbors(lgn_cells, reduction = "pca", dims = 1:15)
lgn_cells <- FindClusters(lgn_cells, resolution = 0.02)

DimPlot(lgn_cells, reduction = "umap", group.by = "ident", alpha = 0.5) 





# save seurat objects
saveRDS(lgn_cells, file="data/lgn_cells.rds")


# get differential genes and save
# lgn_cells <- PrepSCTFindMarkers(lgn_cells,  assay = "SCT", verbose = FALSE)
lgn_cells.markers <- FindAllMarkers(lgn_cells,
                                    test.use = "wilcox",  # Use Wilcoxon rank sum test (you can choose other tests)
                                    min.pct = 0.1,  # Minimum percent of cells expressing the gene
                                    logfc.threshold = 0.25, only.pos = TRUE)  # Log fold change threshold



# Save the markers to a CSV file
write.csv(lgn_cells.markers, file = "data/lgn_cells_markers.csv", row.names = FALSE)


