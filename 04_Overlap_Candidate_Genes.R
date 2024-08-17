setwd("/data/home/jiangminghui/NSCLC")
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(harmony)
library(readr) # 用于读取txt文件
#install.packages("scCustomize")
library(scCustomize)
library(tidyr)
library(gplots)
library(reshape2)
library(GSVA)
library(igraph)
library(ggpubr)
library(RColorBrewer)
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(paletteer)
library(SingleR)
library(celldex)
#install.packages("UpSetR")
library(UpSetR)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(VennDiagram)
#install.packages("venn")
library(ggvenn)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(devtools)
#install_github('ievaKer/aPEAR')
library(DOSE)

library(aPEAR)

#install.packages("ggvenn") # install via CRAN

load('./result/data/WGCNA_keyGene.RData')

load('./result/data/01_significant_markers.RData')

DEG2 <- readRDS(file = './result/data/02_DEG2.rds')

load('./result/data/WGCNA_keyGene.RData')

#load('./result/data/WGCNA_keyGene2.RData')

#load('./result/data/WGCNA_keyGene3.RData')

significant_markers <- significant_markers[significant_markers$pct.1 - significant_markers$pct.2 >0,]

scDEG <- significant_markers %>% rownames()


DEG2 <- subset(DEG2,Sig == "Up") %>% rownames()

length(scDEG)

DEG2<- rownames(DEG2)

length(DEG2)

intersect(intersect(WGCNA_keyGene,scDEG),DEG2) 


venn <- list('DEGs1' =scDEG,
             'MP-MGs' = c(WGCNA_keyGene),
             'DEG2' =  DEG2
               )

str(venn)

pdf('./result/figure/04_Fig4A_Venn_Plot.pdf',width = 6,height = 6)
ggvenn(venn)    
dev.off()


Overlap_genes <- intersect(intersect(venn$DEGs1,venn$`MP-MGs`),venn$DEG2)

Overlap_genes

saveRDS(Overlap_genes,file = './result/data/Overlap_genes.rds')

############################################################
#
#--------------------富集分析————————————————————————————
#
############################################################

Overlap_genes_ENTREZID = bitr(Overlap_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库


head(Overlap_genes_ENTREZID)

ego_ALL <- enrichGO(gene = Overlap_genes_ENTREZID$ENTREZID,
                    OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 1,
                    readable = TRUE) 

pdf('./result/figure/04_Fig4B_GOEnrichment_Plot.pdf',height = 9,width = 6)
dotplot(ego_ALL, showCategory = 10, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale="free")
dev.off()


kk <-  enrichKEGG(
  gene          = Overlap_genes_ENTREZID$ENTREZID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 1,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 1
)
## Reading KEGG annotation online:
##
## Reading KEGG annotation online:

str(kk)
head(kk)
kk<-as.data.frame(kk)
gene_symbols <- bitr(kk$geneID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
pdf('./result/figure/04_Fig4B_KEGGEnrichment_Plot.pdf',height = 8,width = 6)
dotplot(kk, showCategory = 15)
dev.off()


enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = 'CC')

write_csv(as.data.frame(ego_ALL),'./result/data/ego_ALL.csv')
write_csv(as.data.frame(kk),'./result/data/kk.csv')

############################################################
#
#--------------------Goplot富集分析可视化———————————————————
#
############################################################


res_go <- ego_ALL %>% 
  group_by(ONTOLOGY) %>% 
  slice(1:30)

dim(res_go)
res_go <- res_go@result
names(res_go)[c(1,3,7,9)] <- c("Category","Term","adj_pval","Genes")
res_go$Genes <- gsub("\\/",", ",res_go$Genes)


diff_limma <- DEG2[Overlap_genes,]

head(diff_limma)

colnames(diff_limma)

diff_limma$ID <- rownames(diff_limma)

names(diff_limma)[1] <- 'logFC'

go_list <- list(resGo = res_go,
                degLimma = diff_limma
)

circ <- circle_dat(go_list$resGo,go_list$degLimma)

dim(circ)
GOBar(subset(circ, category == 'BP'))
GOBubble(circ, labels = 1)




GOBubble(circ, 
         title = 'Bubble plot', 
         colour = c('orange', 'darkred', 'gold'), 
         display = 'multiple', 
         bg.col = T, #是否显示背景色
         labels = 2) #显示标签的阈值


go_list$genes <- diff_limma[,c("ID","logFC")]
go_list$process <- c("epidermis development",
                     "intermediate filament",
                     "epithelial cell proliferation",
                     "skin development",
                     "keratinocyte differentiation",
                     "perisynaptic extracellular matrix",
                     "npBAF complex","GBAF complex","focal adhesion",
                     "cell-substrate junction","vascular endothelial growth factor receptor binding",
                     "endopeptidase complex"
                     )
length(go_list$process)
chord <- chord_dat(circ, go_list$genes)

d_palettes<- palettes_d_names

col<-c("#3FB8AFFF", "#7FC7AFFF", "#DAD8A7FF", "#FF9E9DFF", "#FF3D7FFF",'#964754FF',
       '#D8AEDDFF', "#BF9BDDFF", "#CB74ADFF","#E69E9CFF", "#FFC3A3FF", "#FBE4C6FF")
mycol<-paletteer_d( "ggsci::default_igv",n=51)

num_rows <- nrow(chord)

pdf('./result/figure/04_Fig4B_GO_Circ.pdf',width = 10,height = 10)
GOChord(chord, space = 0.02, 
        gene.order = 'logFC', 
        gene.space = 0.25, #基因名和图形的距离
        gene.size = 3, #基因名大小
        nlfc = 1, # logfc列的数量？
        lfc.col=c("red","white","blue"), # logfc颜色
        ribbon.col=col,
        lfc.min=-2, #对logfc标准化的最小值
        lfc.max=2, #对logfc标准化的最大值
        #ribbon.col=c("#386CB0","#BEBADA","#D95F02","#8DD3C7","#8DA0CB"),#条目颜色
        border.size=0.8, #条带宽度
        process.label=8, #条目图例标签大小
        limit=c(0,0) 
)

dev.off()

############################################################
#
#--------------------Goplot富集分析可视化（KEGG）——————————
#
############################################################





