###BRCA project###

#@!
#Pipeline& Data handling
library(Seurat)
library(patchwork)
library(dplyr)
library(tidyverse)
library(ggplot2)

#Removing Doublets
library(DoubletFinder)

#DE Analysis
library(DESeq2)

#annotation
library(pheatmap)
library(celldex)
library(SingleR)


#@!
set.seed(1234)
setwd("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test")

#@!
# Load the existing genes.tsv file
genes <- read.delim("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Data\\GSE176078_RAW\\Extractd\\CID4066\\genes.tsv", header = FALSE)
# Create a new data frame with gene IDs (same as names) and gene names
genes_df <- data.frame(V1 = genes$V1, V2 = genes$V1)
# Save the new genes.tsv file
write.table(genes_df, "D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Data\\GSE176078_RAW\\Extractd\\CID4066\\genes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Load the BRCA dataset
BRCA.data <- Read10X(data.dir = "D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Data\\GSE176078_RAW\\Extractd\\CID4066", unique.features = TRUE)
# Initialize the Seurat object with the raw (non-normalized data).
BRCA <- CreateSeuratObject(counts = BRCA.data, project = "GSE176078", min.cells = 3, min.features = 200)
BRCA


#Create a sample column
BRCA$sample<- rownames(BRCA@meta.data)
#Split sample column
BRCA@meta.data<-separate(BRCA@meta.data, col='sample', into= c('patient', 'barcode'), sep= '_')


#@!
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
BRCA[["percent.mt"]] <- PercentageFeatureSet(BRCA, pattern = "^MT-")
length(colnames(BRCA))


# Visualize QC metrics as a violin plot
png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\violin.png", width = 800, height = 600)
VlnPlot(BRCA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(BRCA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BRCA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\scatter.png", width = 800, height = 600)
plot1 + plot2
dev.off()


#Quality Control

#Before QC
length(colnames(BRCA))


BRCA <- subset(BRCA, subset = nFeature_RNA >=200&
#                 nCount_RNA>=1000&
                 nFeature_RNA <= 8000 &
                 percent.mt < 20)

#After QC
length(colnames(BRCA))

#@!
BRCA <- NormalizeData(BRCA)

#@!
BRCA <- FindVariableFeatures(BRCA, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BRCA), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(BRCA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\Variable feature plot.png", width = 1200, height = 600)
plot1 + plot2
dev.off()


#@!
all.genes <- rownames(BRCA)
BRCA <- ScaleData(BRCA, features = all.genes)

#@!
BRCA <- RunPCA(BRCA, features = VariableFeatures(object = BRCA))

# Examine and visualize PCA results a few different ways
print(BRCA[["pca"]], dims = 1:5, nfeatures = 5)

pdf("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\PCA plots.pdf")


VizDimLoadings(BRCA, dims = 1:2, reduction = "pca")
DimPlot(BRCA, reduction = "pca") + NoLegend()
DimHeatmap(BRCA, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(BRCA, dims = 1:15, cells = 500, balanced = TRUE)


dev.off()

png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\Elbowplot.png", width = 800, height = 600)


ElbowPlot(BRCA)


dev.off()


#@1
##Detecting Doublets
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_breast <- paramSweep(BRCA, PCs = 1:10, sct = FALSE)
sweep.stats_breast <- summarizeSweep(sweep.res.list_breast, GT = FALSE)
bcmvn_breast <- find.pK(sweep.stats_breast)


png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\Doublet values.png", width = 800, height = 600)
ggplot(bcmvn_breast, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()
dev.off()

#Setting pK value
pK <- bcmvn_breast %>% 
  filter(BCmetric == max(BCmetric)) %>% 
  select(pK)

pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations<- BRCA@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- BRCA@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(BRCA@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
BRCA <- doubletFinder(BRCA, PCs = 1:10, pN = 0.25, pK = 0.19, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
BRCA@meta.data


png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\Doublet cells.png", width = 800, height = 600)
DimPlot(BRCA, reduction = "pca", group.by = "DF.classifications_0.25_0.17_398")
dev.off()


#Subset doublets
BRCA <- subset(BRCA, subset = DF.classifications_0.25_0.17_398 == 'Singlet')
length(colnames(BRCA))


#@!2
BRCA_neighbors <- FindNeighbors(BRCA, dims = 1:15)
BRCA <- FindClusters(BRCA_neighbors, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(BRCA), 5)

# Look at cluster IDs of the first 5 cells
head(Idents(BRCA), 5)

BRCA <- RunUMAP(BRCA, dims = 1:15)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\UMAP.png", width = 800, height = 600)
DimPlot(BRCA, reduction = "umap")
dev.off()



saveRDS(BRCA, file = "D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\GSE176078.rds") 





BRCA<-readRDS("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\GSE176078.rds")
#Finding Markers for each clusters(((You have to read Vignette )))

#BRCA_markers<- FindAllMarkers(object = BRCA,
#                              logfc.threshold = 0.25,
#                              min.pct = 0.1,
#                             only.pos = FALSE,
#                              test.use = "DESeq2",
#                              slot = "Counts")



BRCA.markers <- FindAllMarkers(BRCA, only.pos = TRUE)
BRCA.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25)


#Top markers for each cluster
BRCA.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(BRCA, features = top10$gene) + NoLegend()
#Feature scatter to see each marker for every cluster


