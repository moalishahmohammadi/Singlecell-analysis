#Cell Trajectory and pseudotime analysis

#Loading required packages
library(ggplot2)
library(tidyverse)
library(Seurat)
library(monocle3)



# Setting directory
setwd("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test")



#Loading data


# ChatGPT recommendations
# Assuming my Seurat object is named 'BRCA'

expression_matrix <- as.matrix(GetAssayData(BRCA, assay = "RNA", layer ="data"))


cell_metadata <- BRCA@meta.data

gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)




#Storing data
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)



#Removing batch effect with cell alignment
cds <- align_cds(cds, alignment_group = "batch")


#Ordering cells
##  Learn a graph
cds <- learn_graph(cds)

##  Order cells
cds <- order_cells(cds)

plot_cells(cds)