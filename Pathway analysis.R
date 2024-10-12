#Pathway analysis
#-------------------------------
#Pathway centric analysis using GSDensity

#Loading Library  
library(gsdensity)
library(ggplot2) # for plotting
library(reshape2)
library(msigdbr) # for gathering gene sets
library(Seurat)
library(SeuratData)
library(future) # for parallel computing
# library(future.apply) # for parallel computing

# use GO_BP gene sets 
# Conver the format to 'list'

mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
mdb_c5_bp <- mdb_c5[mdb_c5$gs_subcat == "GO:BP", ]
gene.set.list <- list()
for (gene.set.name in unique(mdb_c5_bp$gs_name)){
  gene.set.list[[gene.set.name]] <- mdb_c5_bp[mdb_c5_bp$gs_name %in% gene.set.name, ]$gene_symbol
}


#Loading data

BRCA<- readRDS("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\annotated.rds")

# from SeuratData
data("BRCA")
# Normalize the data
BRCA <- NormalizeData(BRCA)

# compute cell/gene embeddings
ce <- compute.mca(object = BRCA)
# find gene sets with differential density 
res <- compute.kld(coembed = ce, 
                   genes.use = intersect(rownames(ce), rownames(BRCA)), # this intersection is to select only genes, not cells. 
                   n.grids = 100, 
                   gene.set.list = gene.set.list,
                   gene.set.cutoff = 3,
                   n.times = 100)                   
gene.set.deviated <- res[res$p.adj < 0.05, ]$gene.set

# build nearest neighbor graph
cells <- colnames(BRCA)
el <- compute.nn.edges(coembed = ce, nn.use = 300)

# get label propagation probability for each cell of a gene set 'GOBP_B_CELL_ACTIVATION'
cv <- run.rwr(el = el, gene_set = gene.set.list[["GOBP_B_CELL_ACTIVATION"]], cells = cells)
# get a 'positive' or 'negative' label for each cell of this gene set
cl <- compute.cell.label(cv)

# Compute the UMAP coordinates for visualization
BRCA <- BRCA %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:20)
# plot the cell annotation
DimPlot(BRCA,
        group.by = "seurat_annotations",
        label = T,
        raster = T)