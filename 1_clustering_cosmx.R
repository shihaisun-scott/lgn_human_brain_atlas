

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
lgn_cells_filtered <- subset(lgn_cells, features = keepgen)


# lgn_cells <- NormalizeData(object = lgn_cells)
lgn_cells <- FindVariableFeatures(object = lgn_cells, selection.method = "vst", nfeatures = 3000)
lgn_cells <- ScaleData(object = lgn_cells, verbose = FALSE)
lgn_cells <- RunPCA(object = lgn_cells, npcs = 30, verbose = FALSE)
ElbowPlot(lgn_cells, ndims = 30)

lgn_cells <- RunUMAP(object = lgn_cells, reduction = "pca", dims = 1:20)
lgn_cells <- FindNeighbors(lgn_cells, reduction = "pca", dims = 1:20)
lgn_cells <- FindClusters(lgn_cells, resolution = 0.2)

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


