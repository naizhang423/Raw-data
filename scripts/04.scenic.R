# https://cloud.tencent.com/developer/article/1692240
# https://zhuanlan.zhihu.com/p/358986392
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('/home/pub252/projects/codes/mg_base.R')

library(Seurat)
# library(SeuratObject)
library(Matrix)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(Seurat)
library(doParallel)
library(SCopeLoomR)
setwd
mydata = readRDS
set.seed(1234)
subcell <- sample(colnames(mydata),1000)
scRNAsub <- mydata[,subcell]

exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
cellInfo = mydata@meta.data[, c(11, 3, 2)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')

### Initialize settings
library(SCENIC)
# cisTarget_databases 
# 
mydbDIR <- ""
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", dbDir=mydbDIR, dbs=mydbs, nCores=4, datasetTitle="T cells") 
saveRDS(scenicOptions, "int/scenicOptions.rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene=3*0.01*ncol(exprMat), minSamples=ncol(exprMat)*0.01)
exprMat_filtered <- exprMat[genesKept, ]
## 
runCorrelation(exprMat_filtered, scenicOptions)
## 
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts=20)
## 
### Build and score the GRN
exprMat_log <- log2(exprMat+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top10perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")        # choose settings

####################################################################################################################
# 
# 
##
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
# mydata <- AddMetaData(mydata, AUCmatrix) #
scRNAsub <- AddMetaData(scRNAsub, AUCmatrix) #
##################################################################
## 

## 

## 

# df = read.table("./TF_AUC_heatmap.txt", header=T, row.names=1, sep="\t")
# pheatmap(df, scale="row", angle_col="45", color=colorRampPalette(c("#0da9ce", "white", "#e74a32"))(500), clustering_method="median", border_color="white")


head(mydata@meta.data)
df=data.frame(group=scRNAsub@meta.data$Samples2,scRNAsub@meta.data[,colnames(AUCmatrix)])

head(df)
pvalue <- c()
for (i in 2:ncol(df)) {
  df_PT <- df[df$group=='PT',i]
  df_Li <- df[df$group=='Li',i]
  p <- wilcox.test(df_PT,df_Li)$p.value
  pvalue <- c(pvalue,p)
}
pvalue <- c(0,pvalue)
df <- rbind(df,pvalue)
df <- df[,which(df[1001,]<0.001)]
df <- df[-1001,]
df=melt(df)

ggplot(df, aes(x=variable, y=value, fill=group))+
  scale_fill_manual(values=c("dodgerblue", "hotpink"))+
  geom_boxplot(outlier.size=0.1, width=0.3)+
  theme_bw()+
  stat_compare_means(aes(group=group), label="p.signif", method="wilcox")+
  theme(axis.text.x=element_text(angle=90, hjust=1, face="bold", size=10), axis.text.y=element_text(face="bold", size=10))
ggsave("boxplot.pdf",width = 10,height = 6)

##FOS_extended_644g,JUNB_extended_1033g


diff_pathway<-function(dat,group){
  dat=data.frame(cluster=group,t(dat))
  #gr=c('High','Low')
  gr=names(table(group))
  dat1=dat[dat$cluster==gr[1],-1]
  dat2=dat[dat$cluster==gr[2],-1]
  # dat3=dat[dat$cluster==gr[3],-1]
  pathway=unique(colnames(dat)[-1])
  p_vale=data.frame()
  for (i in pathway){
    # x=c(dat1[,i],dat2[,i],dat3[,i])
    # g= factor(rep(names(table(group)), c(nrow(dat1), nrow(dat2), nrow(dat3))),
    #           labels = names(table(group)))
    # dd1=kruskal.test(x,g)$p.value
    dd1=wilcox.test(dat1[,i],dat2[,i])$p.value
    p_vale=rbind(p_vale,data.frame(pathway=i,p.value=dd1))
  }
  return(p_vale)
}
aucell.diff<-diff_pathway(dat=t(scRNAsub@meta.data[,c("FOS_extended_644g","JUNB_extended_1033g")]),group=scRNAsub@meta.data$Samples2)

tf_regulon=readRDS("int/3.1_regulons_forAUCell.Rds")
FOS_enrich_res=mg_clusterProfiler(genes = tf_regulon$`FOS_extended (644g)`)
enrichplot::dotplot(FOS_enrich_res$KEGG)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
ggsave("FOS_enrich_KEGG.pdf",height = 6,width = 8)
enrichplot::dotplot(FOS_enrich_res$GO_BP)+ggtitle('Biological Process')+
  theme(text=element_text(family = 'Times'))
ggsave("FOS_enrich_GOBP.pdf",height = 6,width = 8)
write.csv(FOS_enrich_res$GO_BP,'FOS_enrich_GOBP.csv')

JUNB_enrich_res=mg_clusterProfiler(genes = tf_regulon$`JUNB_extended (1033g)`,keytype = 'SYMBOL')
enrichplot::dotplot(JUNB_enrich_res$KEGG)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
ggsave("JUNB_enrich_KEGG.pdf",height = 6,width = 8)
enrichplot::dotplot(JUNB_enrich_res$GO_BP)+ggtitle('Biological Process')+
  theme(text=element_text(family = 'Times'))
ggsave("JUNB_enrich_GOBP.pdf",height = 6,width = 8)
write.csv(JUNB_enrich_res$GO_BP,'JUNB_enrich_GOBP.csv')
