
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
library(openxlsx)
library(gridExtra)

# load differential markers
lgn_cells.markers <- read.csv("data/lgn_cells_markers.csv")

# Load the processed Seurat object
lgn_cells <- readRDS("data/lgn_cells.rds")
plot_list = list()

# plot UMAP
umap_plot <- DimPlot(lgn_cells, reduction = "umap", group.by = "ident") + 
  ggtitle("") + 
  theme(aspect.ratio = 1,
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
plot_list[[1]] <- umap_plot


# get gene names
gene_names <- lgn_cells@assays$RNA@meta.features$feature_name
# gene_names <- lgn_cells@assays$RNA@meta.features$Gene # this produces the same results

# plot dot plots for specific LGN markers
# specific_lgn_genes <- c(
#   "CALB1", "PVALB", "VGLUT", "NEFH", "TFAP2B", "RAX", "BDNF", "GAD1", "SYP", "CAMK2A"
# )

specific_lgn_genes <- c(
  "CALB1", "PVALB", "FOXP2", "SLC17A7", "TFAP2B", "GAD1", "CAMK2A",
  "NPSR1", "RGS4", "NRCAM", "SLC17A6", "NTRK2", "PPP1R1B", "STXBP6", "CACNA1C", "GRM2", "KCNQ3",  # M cells
  "GRIN2B", "CAMK2A", "HOMER1", "SYN1", "DLGAP1", "SLC6A1", "GAD1", "SHANK2", "GABRA1", "SV2A",  # P cells
  "CALB1", "SLC17A7", "NTRK3", "MEF2C", "LAMP5", "PVALB", "GRIK1", "NPY", "VIP", "STMN2"  # K cells
)

specific_lgn_genes <- unique(specific_lgn_genes)

specific_genes_ids <- rownames(lgn_cells)[gene_names %in% specific_lgn_genes]

# Plot dot plot with specific gene IDs
plot_list[[2]] <-DotPlot(lgn_cells, features = specific_genes_ids, group.by = "seurat_clusters", scale.min = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels = specific_lgn_genes)

# do the same as above but match with annotation file
gtf_file <- "data/gencode.v44.primary_assembly.annotation.gtf"
gtf_data <- import(gtf_file, format = "gtf")
specific_genes_ids2 = c()
n_unique = c()
new_gene_names <- c()
for (ii in specific_lgn_genes) {
  unique_geneids <- length(unique(gtf_data$gene_id[grep(ii, gtf_data$gene_name)]))
  if (unique_geneids == 1){
    specific_genes_ids2 <- c(specific_genes_ids2, unique(gtf_data$gene_id[grep(ii, gtf_data$gene_name)]))
    n_unique <- c(n_unique, length(unique(gtf_data$gene_id[grep(ii, gtf_data$gene_name)])))
    new_gene_names <- c(new_gene_names, ii)
  }
}

specific_genes_ids2_clean <- sub("\\..*$", "", specific_genes_ids2)

plot_list[[2]] <- DotPlot(lgn_cells, features = specific_genes_ids2_clean, group.by = "seurat_clusters", scale.min = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels = new_gene_names)



# save figures into a pdf
pdfname <- paste0("analysis_output/lgn_cells.pdf")
# pdf(pdfname, height = 12, width = 18)
# grid.arrange(umap_plot, heat_map, dot_plot, dot_plot2, nrow = 2, ncol = 2)
# dev.off()

pdf(pdfname, height = 7, width = 11)
for (plot in plot_list) {
  print(plot)
}
dev.off()

