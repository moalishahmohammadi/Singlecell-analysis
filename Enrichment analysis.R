#Enrichment analysis using cluster profiler


#Required packages for loading data
library(Seurat)
library(dplyr)

#Loading Seurat object
seurat_obj<- readRDS("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\annotated.rds")


# Extract gene expression data ---- Can i do it for a specific cluster or cell subset?
cluster_data <- seurat_obj[["RNA"]]$counts

# Convert to a data frame
cluster_df <- data.frame(gene = rownames(cluster_data), expression = rowMeans(cluster_data))

# Sort by expression
cluster_df <- cluster_df %>% arrange(desc(expression))

#Performing cluster profiler
library(clusterProfiler)
data(geneList, package="DOSE")
## fold change > 2 as DE genes
de <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
## use simplify to remove redundant terms
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)


#GSEA using KEGG pathways and modules
kk <- gseKEGG(geneList, organism = "hsa") 














