#Cell-Cell communication

# loading data
BRCA<- readRDS("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\annotated.rds")

#Loading libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T) 


data.input<-BRCA[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `BRCA[["RNA"]]$data`
labels <- Idents(BRCA)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels


# Creating the cell chat object

cellchat <- createCellChat(object = BRCA, group.by = "ident", assay = "RNA")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


#Setting the database for ligand-receptor interaction

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)


# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling


# set the used database in the object
cellchat@DB <- CellChatDB.use



#Preprocessing the expression data 


# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692
# Start measuring execution time

ptm <- Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 13.20763
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)


#Computing communication probability

cellchat <- computeCommunProb(cellchat, type = "triMean")
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-02-14 00:32:35.767285]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-02-14 00:33:13.121225]"


#Filtering number of cells for each group

cellchat <- filterCommunication(cellchat, min.cells = 10)


#Extracting cellular communication etwork as the form of dataframe

df.net <- subsetCommunication(cellchat, slot.name = "netP")


#Extracting cellular communication at signalling pathways level
cellchat <- computeCommunProbPathway(cellchat)


#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 38.73308


#Visualizing interactions
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\cell-cell.png", width = 2000, height = 1600)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\cell-cell2.png", width = 2000, height = 1600)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector.
png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\cell-cell3.png", width = 2000, height = 1600)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()
# Circle plot
par(mfrow=c(1,1))
png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\cell-cell4.png", width = 2000, height = 1600)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()



#Khafane

# Chord diagram
par(mfrow=c(1,1))
png("D:\\Bioinformatics\\Projects\\CAF's\\One by one\\Test\\cell-cell5.png", width = 2000, height = 1600)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

