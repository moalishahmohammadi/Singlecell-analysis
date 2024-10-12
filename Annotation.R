
#Require libraries
library(pheatmap)
library(celldex)
library(SingleR)


setwd("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test")
#Annotation
#Manual
png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\Feature scatter for markers.png", width = 800, height = 600)

FeaturePlot(BRCA, features = c("CCL5", "IL7R", "SFRP2", "NDUFA4L2",
                               "RAMP2", "KRT19", "CTLA4", "CXCL12",
                               "SLPI", "KRT14", "STC2", "RRM2"))

dev.off()


new.cluster.ids <- c("Mesenchymal cell", "CD4+ T cell", "Luminal progenitor cell", "Perivascular cell",
                     "Cancer cell", "Fibroblast", "Regulatory T(Treg) cell","Macrophage",
                     "Fibroblast", "Luminal progenitor cell", "Basal epithelial cell","Mature luminal cell", "Proliferating T cell
")
names(new.cluster.ids) <- levels(BRCA)
BRCA <- RenameIdents(BRCA, new.cluster.ids)
DimPlot(BRCA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


library(ggplot2)
plot <- DimPlot(BRCA, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test.png", height = 7, width = 12, plot = plot)







#Automatic
##Annotation USING SingleR
# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the test dataset

BRCA_counts <- GetAssayData(BRCA, layer ='counts')

pred <- SingleR(test = BRCA_counts,
                ref = ref,
                labels = ref$label.main)

BRCA$singleR.labels <- pred$labels[match(rownames(BRCA@meta.data), rownames(pred))]
DimPlot(BRCA, reduction = 'umap', group.by = 'singleR.labels')

# Annotation diagnostics ----------
# ...Based on the scores within cells -----------
pred
pred$scores

plotScoreHeatmap(pred)

# ...Based on deltas across cells ----------

plotDeltaDistribution(pred)

# ...Comparing to unsupervised clustering ------------

tab <- table(Assigned=pred$labels, Clusters=BRCA$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))




#Save progress

saveRDS(BRCA, "D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\annotated.rds")
