# 
# 
# 
devtools::install_github("junjunlab/GseaVis")
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GseaVis)
library(Seurat)
setwd
# 1
mydata = readRDS
mydata = subset(mydata, cell_type=="Cancer stem cell")
mydata$Samples2 <- substr(mydata$Samples,1,2)
data2 = FindMarkers(mydata, group.by="Samples2", ident.1="PT", ident.2="Li", logfc.threshold=0.25, min.pct=0)
write.table(data2, "./TNK_treatment_control_DEG0.25.txt", col.names=T, row.names=T, quote=F, sep="\t")

# symbol entrez ID：
symbol <- rownames(data2)
entrez <- bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(entrez)
#    SYMBOL ENTREZID
# 11   A1BG        1
# 12   A1CF    29974
# 13    A2M        2
# 14 A4GALT    53947
# 15  A4GNT    51146
# 16   AAAS     8086

# (entrez+log2FC)：
genelist <- data2$avg_log2FC
names(genelist) <- rownames(data2)
#
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist), entrez[,1]), 2]
genelist = sort(genelist, decreasing=T)

##############################################################################################
# 
R.utils::setOption("clusterProfiler.download.method", "auto") ##
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",
  # minGSSize = 10,
  # maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
)
KEGG_ges2 <- setReadable(KEGG_ges,
                         OrgDb = org.Hs.eg.db,
                         keyType="ENTREZID")
#
KEGG_ges@result$core_enrichment[1]
#
KEGG_ges2@result$core_enrichment[1]
save(KEGG_ges,KEGG_ges2,file = c('GSEA_KEGG.Rdata'))
write.table(KEGG_ges2@result,"KEGG_result.txt",quote = F,sep = "\t")

R.utils::setOption("clusterProfiler.download.method", "auto") ##
GO_ges <- gseGO(
  geneList=genelist,
  ont="BP",
  OrgDb=org.Hs.eg.db,
  keyType="ENTREZID",
  # minGSSize=10,
  # maxGSSize=500,
  pvalueCutoff=1,
  pAdjustMethod="BH",
  verbose=FALSE,
)
GO_ges2 <- setReadable(GO_ges,
                       OrgDb = org.Hs.eg.db,
                       keyType="ENTREZID")
#
GO_ges@result$core_enrichment[1]
#
GO_ges2@result$core_enrichment[1]
save(GO_ges,GO_ges2,file = c('GSEA_GO.Rdata'))
write.table(GO_ges2@result, "./GOBP_result.txt", col.names=T, row.names=T, quote=F, sep="\t")
#################################################################################
# 
########----------------------------------------------------------------------####
#
BP <- KEGG_ges2@result
BP <- na.omit(BP)
BP<- BP[order(BP$NES,decreasing = F),]
BP$Description<-factor(BP$Description,levels=unique(BP$Description))

p1<-ggplot(data = BP, aes(x=NES, y=Description, fill=pvalue))+
  geom_bar(stat="identity")+ 
  scale_fill_gradient(low="#99CC99",high ="grey80")+
  theme_minimal()

p1
ggsave('kegg_bar.pdf',height = 7,width = 7)
#2.
result <- KEGG_ges2@result
View(result) #
#
##
which(result$Description=="Protein export")
core <- KEGG_ges2@result$core_enrichment[87]
length(core)#
core
#
core_genes <- str_split(core ,'/')[[1]]
core_genes
gseaNb(
  object = KEGG_ges2,
  geneSetID = KEGG_ges2@result$ID[87],
  addPval = T,
  pvalX = 0.95,
  pvalY = 0.8,
  newGsea = T,
  addPoint = F,
  newCurveCol = c("#001871","#b9b4ad", "#f99f1c"),
  newHtCol = c("#001871", "white", "#f99f1c"),
  addGene = core_genes, #
  kegg = T,
  geneCol = '#4d4d4d'
)
ggsave("Protein export.pdf",height = 6,width = 6)
##########------------------------------------------------

DotPlot(mydata, group.by ="Samples2"  ,features = core_genes) + 
  scale_colour_gradient2(low = '#0099CC', mid = "lightgrey", high ="#FF9900" )+
  scale_y_discrete(labels = c("PT", "Li"))+
  labs(title='Protein export')+
  coord_flip()+theme_bw()+theme(axis.text.x = element_text(size=10, angle=0, hjust=0.5),text = element_text(family = "serif",face = "bold",color = "black"))
ggsave("dot_gene_Pentose phosphate pathway.pdf",width=4,height =6)


##########################################################################
# 
library(ggplot2)
df = read.table("./GSEA_GOBP_result_barplot.txt", header=T, sep="\t")
df <- result[result$pvalue<0.05,]
ggplot(df, aes(x=NES, y=reorder(Description, NES), color=pvalue))+
  geom_point(aes(size=abs(NES)))+
  geom_segment(aes(x=0, xend=NES, y=Description, yend=Description), color="black", linewidth=1)+
  scale_color_gradient(low="hotpink1", high="steelblue1")+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=12), axis.text.y=element_text(face="bold", size=8), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=12))
ggsave("Lollipop diagram.pdf",width=8,height =6)




