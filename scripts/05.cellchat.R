rm(list = ls())
# 
# devtools::install_github("sqjin/CellChat")
library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(svglite)
options(stringsAsFactors=FALSE)
setwd
#####  
mydata = readRDS
mydata_PT <- subset(mydata,Samples2=="PT")
mydata_Li <- subset(mydata,Samples2=='Li')
cellchat = createCellChat(object=mydata, group.by="cell_type", meta=mydata@meta.data)
groupSize = as.numeric(table(cellchat@idents)) #

#####  
CellChatDB = CellChatDB.human     #CellChatDB = CellChatDB.mouse

#####  
# > unique(CellChatDB$interaction$annotation)
# [1] "Secreted Signaling"	"ECM-Receptor"	"Cell-Cell Contact"  #

CellChatDB.use <- subsetDB(CellChatDB, search="Secreted Signaling")
cellchat@DB <- CellChatDB.use

#####  
cellchat = subsetData(cellchat)
future::plan("multisession", workers=4)
## 
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)  #

## 
cellchat = projectData(cellchat, PPI.human)
## 

#####################################################################################
### 
cellchat = computeCommunProb(cellchat, raw.use=FALSE, population.size=TRUE)
# 
cellchat = filterCommunication(cellchat, min.cells=10)
df.net = subsetCommunication(cellchat)

#####  
cellchat = computeCommunProbPathway(cellchat)
df.netp = subsetCommunication(cellchat, slot.name="netP")

#####################################################################################
#######  
## 
cellchat = aggregateNet(cellchat)
## 
groupSize = as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight=groupSize, weight.scale=T, label.edge=T, title.name="Number of interactions", arrow.size=1, color.use=c("hotpink", "deepskyblue", "forestgreen"), edge.label.cex=1.2)
netVisual_circle(cellchat@net$weight, vertex.weight=groupSize, weight.scale=T, label.edge=T, title.name="Number of strength", arrow.size=1, color.use=c("hotpink", "deepskyblue", "forestgreen"))







########
cellchat_PT = createCellChat(object=mydata_PT, group.by="cell_type", meta=mydata_PT@meta.data)
cellchat_Li = createCellChat(object=mydata_Li, group.by="cell_type", meta=mydata_Li@meta.data)

cellchat_PT@DB = CellChatDB.human
cellchat_PT = subsetData(cellchat_PT)
cellchat_PT = identifyOverExpressedGenes(cellchat_PT)
cellchat_PT = identifyOverExpressedInteractions(cellchat_PT)
cellchat_PT = projectData(cellchat_PT, PPI.human)
cellchat_PT = computeCommunProb(cellchat_PT, raw.use=FALSE, population.size=TRUE)
cellchat_PT = filterCommunication(cellchat_PT, min.cells=10)
cellchat_PT = computeCommunProbPathway(cellchat_PT)
cellchat_PT = aggregateNet(cellchat_PT)

cellchat_Li@DB = CellChatDB.human
cellchat_Li = subsetData(cellchat_Li)
cellchat_Li = identifyOverExpressedGenes(cellchat_Li)
cellchat_Li = identifyOverExpressedInteractions(cellchat_Li)
cellchat_Li = projectData(cellchat_Li, PPI.human)
cellchat_Li = computeCommunProb(cellchat_Li, raw.use=FALSE, population.size=TRUE)
cellchat_Li = filterCommunication(cellchat_Li, min.cells=10)
cellchat_Li = computeCommunProbPathway(cellchat_Li)
cellchat_Li = aggregateNet(cellchat_Li)

# 
cellchat_list = list(PT=cellchat_PT, Li=cellchat_Li)
cellchat_combine = mergeCellChat(cellchat_list, add.names=names(cellchat_list), cell.prefix=T)
# 
compareInteractions(cellchat_combine, show.legend=T, group=c('PT', 'Li'), measure="count", color.use=c("#ABC567", "#FF99BB"))
ggsave("cellchat_count.pdf",width = 8,height = 6)
# 
rankNet(cellchat_combine, mode="comparison", stacked=T, do.stat=T, color.use=c("#ABC567", "#FF99BB"))
ggsave("cellchat_strength.pdf",width = 8,height = 6)

########
pdf("cellchat_PT.pdf",width = 6,height = 6)
groupSize = as.numeric(table(cellchat_PT@idents))
netVisual_circle(cellchat_PT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("cellchat_Li.pdf",height = 6,width = 6)
groupSize = as.numeric(table(cellchat_Li@idents))
netVisual_circle(cellchat_Li@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("cellchat_Li.pdf",height = 6,width = 6)
groupSize = as.numeric(table(cellchat_Li@idents))
netVisual_circle(cellchat_Li@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()





### 
mat = cellchat_PT@net$count    ### or mat = cellchat@net$weight
par(mfrow=c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)){
  mat2 = matrix(0, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2, vertex.weight=groupSize, weight.scale=T, arrow.width=0.2, arrow.size=0.1, edge.weight.max=max(mat), title.name=rownames(mat)[i])
}

####TF####
## 
## 
netAnalysis_contribution(cellchat_PT, signaling=cellchat_PT@netP$pathways)
cellchat_PT@netP$pathways
## 
pairLR = extractEnrichedLR(cellchat_PT, signaling=cellchat_PT@netP$pathways, geneLR.return=FALSE)
#########################################################################################
## 
levels(cellchat_PT@idents)  #
## 
netVisual_bubble(cellchat_PT,remove.isolate=FALSE)
ggsave("TF_PT.pdf",height = 10,width = 8)
## 
# netVisual_bubble(cellchat, sources.use=c(3,5,7,8,9), targets.use=c(1,2,4,6), signaling=c("CCL", "CXCL"), remove.isolate=FALSE)
# 
# netVisual_bubble(cellchat, signaling=c("IL2"), remove.isolate=FALSE)+theme(axis.text.x=element_text(angle=15, hjust=1))
# groupSize <- as.numeric(table(cellchat@idents))


netAnalysis_contribution(cellchat_Li, signaling=cellchat_Li@netP$pathways)
cellchat_Li@netP$pathways
## 
pairLR = extractEnrichedLR(cellchat_Li, signaling=cellchat_Li@netP$pathways, geneLR.return=FALSE)
#########################################################################################
## 
levels(cellchat_Li@idents)  #
## 
netVisual_bubble(cellchat_Li,remove.isolate=FALSE)
ggsave("TF_Li.pdf",height = 10,width = 8)
## 
# netVisual_bubble(cellchat, sources.use=c(3,5,7,8,9), targets.use=c(1,2,4,6), signaling=c("CCL", "CXCL"), remove.isolate=FALSE)
# 
# netVisual_bubble(cellchat, signaling=c("IL2"), remove.isolate=FALSE)+theme(axis.text.x=element_text(angle=15, hjust=1))
# groupSize <- as.numeric(table(cellchat@idents))

