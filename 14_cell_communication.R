library(CellChat)
library(tidyverse)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

Idents(GSE117570_SC) <- 'Subcluster'
GSE117570_SCe_marke <- FindAllMarkers(GSE117570_SC,only.pos = T,min.pct = 0.25, 
                                      logfc.threshold = 1,
                                      min.diff.pct = 0.25)

top_markers <- GSE117570_SCe_marke %>% group_by(cluster) %>% 
  top_n(n = 3, wt = avg_log2FC) 

DefaultAssay(GSE117570_SC)
pdf('./result/figure/14_Fig14A_DoHeatmap.pdf',width = 8,height = 6)
DoHeatmap(object = GSE117570_SC,group.by ='Subcluster' ,slot = 'count',features = top_markers$gene)
dev.off()
str(GSE117570_SC)

cellchat <- createCellChat(object = GSE117570_SC,
                           meta = GSE117570_SC@meta.data,
                           group.by = "Subcluster")
cellchat
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 

cellchat <- identifyOverExpressedGenes(cellchat)
#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat)

#project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)
#pdf('./result/figure/14B')

pdf('./result/figure/14_fig14B_cell_communication.pdf')
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()

pdf('./result/figure/14_fig14C_cell_communication.pdf')

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()
#指定TGFb和SPP1两个信号通路
#指定受体-配体细胞类型
pdf('./result/figure/14_fig14D_cell_communication.pdf')

netVisual_bubble(cellchat, sources.use = c(3,5), 
                 targets.use = c(1,2,4,6), remove.isolate = FALSE)
dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat,
                                          slot.name = "netP")
# 所有信号通路
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
# 指定信号通路
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MHC-I"))
gg1 + gg2
# 所有信号通路
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
# 指定信号通路
ht3 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-I", "MHC-II"), pattern = "outgoing")
ht4 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-I", "MHC-II"), pattern = "incoming")
ht3 + ht4

selectK(cellchat, pattern = "outgoing")

cellchat <- identifyCommunicationPatterns(cellchat,
                                          k = 2, #根据selectK图形确定具体K值
                                          pattern = "outgoing"
)
pdf('./result/figure/14_fig14E_cell_communication1.pdf')
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

pdf('./result/figure/14_fig14E_cell_communication1.pdf')
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

selectK(cellchat, pattern = "incoming")
cellchat <- identifyCommunicationPatterns(cellchat,
                                          k = 3, #根据selectK图形确定具体K值
                                          pattern = "incoming"
)
pdf('./result/figure/14_fig14E_cell_communication2.pdf')
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
# dot plot

netAnalysis_dot(cellchat, pattern = "incoming")



















