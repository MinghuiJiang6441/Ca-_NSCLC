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
library(DESeq2)
library(edgeR)
library(igraph)
library(tinyarray)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

library(paletteer)

library(sva)


#解压zip格式文件
#utils::unzip('./data/LUAD_merge.mRNA.TPM.txt.zip', exdir = "./data")#unzip
#utils::unzip('./data/LUSC_merge.mRNA.TPM.txt.zip', exdir = "./data")#unzip

LUSC_TPM <- read.delim(file = './data/LUSC_merge.mRNA.TPM.txt', header = TRUE)

LUAD_TPM <- read.delim(file = './data/LUAD_merge.mRNA.TPM.txt', header = TRUE)

# 定义一个函数来处理数据框
process_dataframe <- function(df) {
  # 1. 去掉ENGSID列
  df <- df %>% select(-ENGSID)
  
  # 2. 将Symbol列变成行名（注意不要有重复）
  df <- df %>% distinct(Symbol, .keep_all = TRUE) %>% column_to_rownames(var = "Symbol")
  
  # 3. 列名的“.”变为“-”
  colnames(df) <- gsub("\\.", "-", colnames(df))
  
  # 4. 删除type列
  df <- df %>% select(-type)
  
  return(df)
}


LUAD_TPM_processed <- process_dataframe(LUAD_TPM)

LUSC_TPM_processed <- process_dataframe(LUSC_TPM)


LUSC_TPM_processed[1:5,1:5]
# 找到共有基因
common_genes <- intersect(rownames(LUAD_TPM_processed), rownames(LUSC_TPM_processed))
length(common_genes)


# 提取共有基因的数据
LUAD_common <- LUAD_TPM_processed[common_genes, ]
LUSC_common <- LUSC_TPM_processed[common_genes, ]

# 合并数据框
TCGA_NSCLC_exp <- cbind(LUAD_common, LUSC_common)

range(TCGA_NSCLC_exp)

#TCGA_NSCLC_exp <- log2(TCGA_NSCLC_exp+0.1)

#TCGA_NSCLC_exp <- log2(TCGA_NSCLC_exp+1)

range(TCGA_NSCLC_exp)
#keepData <- rowSums(TCGA_NSCLC_exp>0) >= floor(0.75*ncol(TCGA_NSCLC_exp))

#TCGA_NSCLC_exp <- TCGA_NSCLC_exp%>% mutate(across(everything(), ~ . * 1e6 / sum(.)))#转化为count
#non_integer_values <- TCGA_NSCLC_exp[!apply(TCGA_NSCLC_exp, 2, function(x) all(x == floor(x)))]
#TCGA_NSCLC_exp <- round(TCGA_NSCLC_exp)
############################################################
#
#-------------------------消除批次效应--------————————————————
#
############################################################
#因为是两套数据集合并到一起的 所以消除批次
colData = data.frame(
  condition = as.factor(c(colnames(LUAD_common),colnames(LUSC_common) )),
  type = as.factor(c(rep('LUAD',ncol(LUAD_common)),rep('LUSC',ncol(LUSC_common)))))

colData$condition %in% colnames(TCGA_NSCLC_exp) %>% any #TRUE

rownames(colData) <- colnames(TCGA_NSCLC_exp)

draw_pca(exp = TCGA_NSCLC_exp, group_list = factor(colData$type))#批次效应可视化

expr_combat <- ComBat(dat = TCGA_NSCLC_exp, batch = colData$type)

expr_combat[1:4,1:4]
expr_combat <- as.data.frame(expr_combat)

pdf('./result/figure/02_SupFig_DEG2_Batch_Effect.pdf',width = 6,height = 5)
draw_pca(exp = expr_combat, group_list = factor(colData$type))
dev.off()
dim(expr_combat)

range(expr_combat)

############################################################
#
#-------------------------差异表达基因--------————————————————
#
############################################################

group=ifelse(substring(colnames(expr_combat),14,15) == "11","Normal","Tumor")
expr_combat[expr_combat < 0] <- 0
y <- DGEList(counts=expr_combat,group=group) ##Remove rows consistently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##Perform TMM normalization and transfer to CPM (Counts Per Million)

y <- calcNormFactors(y,method="TMM")

count_norm=cpm(y)

count_norm<-as.data.frame(count_norm)

#Run the Wilcoxon rank-sum test for each gene
pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),group)
  p=wilcox.test(gene~group, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")

#Calculate the fold-change for each gene:
group <- as.factor(group)

conditionsLevel<-levels(group)

dataCon1=count_norm[,c(which(group==conditionsLevel[1]))]

dataCon2=count_norm[,c(which(group==conditionsLevel[2]))]

foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

#Output results based on FDR threshold
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
fdrThres=0.05

#outRst$padj <- p.adjust(outRst$pValues, method = "BH")

head(outRst)

saveRDS(TCGA_NSCLC_exp,file = "./data/TCGA_NSCLC_exp.rds") #save
#没有log 去了批次


#TCGA_NSCLC_exp <- readRDS(file = "./data/TCGA_NSCLC_exp.rds")
############################################################
#
#-----------------------------------火山图------------------
#
############################################################

cut_off_padj =0.05 #设置padj的阈值

cut_off_log2FC =0.5 #设置log2FC的阈值

outRst$Sig = ifelse(outRst$FDR < cut_off_padj &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    abs(outRst$log2foldChange) >= cut_off_log2FC,  #abs绝对值
                  ifelse(outRst$log2foldChange > cut_off_log2FC ,'Up','Down'),'no')

table(outRst$Sig)
head(outRst)

#显示标签的基因
# 按FDR排序并选择前10个
labeldata <- subset(outRst, abs(log2foldChange) > 1 & FDR < 0.01)
labeldata$label <- rownames(labeldata)

# 按FDR排序并在每个Sig组中选择前10个
labeldata <- labeldata[order(labeldata$Sig, -log10(labeldata$FDR),decreasing = TRUE), ]
labeldata <- do.call(rbind, lapply(split(labeldata, labeldata$Sig), function(x) head(x, 10)))

pdf('./result/figure/02_Fig2D_DEG2_Volcano_Plot.pdf',width = 9,height = 8)
ggplot(outRst, aes(x = log2foldChange, y = -log10(FDR), colour = Sig)) + 
  geom_point(alpha = 0.65, size = 2) +  
  scale_color_manual(values = c("#546de5", "#d2dae2", "#ff4757")) + 
  xlim(c(-7, 8)) +  
  geom_vline(xintercept = c(-cut_off_log2FC, cut_off_log2FC), lty = 4, col = "black", lwd = 0.8) + 
  geom_hline(yintercept = cut_off_padj, lty = 4, col = "black", lwd = 0.8) +  
  labs(x = "log2FC", y = "-log10padj") +  
  ggtitle("Volcano Plot") +  # 更改标题为英文
  theme_bw() +
  geom_text_repel(data = labeldata, aes(label = label), max.overlaps = 30, color = "black") 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right", 
        legend.title = element_blank()
  )
dev.off()

DEG2 <- outRst 
saveRDS(DEG2,file = './result/DEG2.rds')

############################################################
#
#------------------------------环形热图------------------
#
############################################################
#exarct DEG
TCGA_NSCLC_exp<-readRDS("./data/TCGA_NSCLC_exp.rds")

DEG2 <- readRDS('./result/DEG2.rds')
#!grepl("11$", substr(colnames(TCGA_NSCLC_exp), 14, 15))%>% table

#use.gene <- subset(x = DEG2, Sig  != 'no') %>% rownames()
use.gene_df <- subset(DEG2,Sig=='Up')

use.gene_df <- use.gene_df[order(use.gene_df$Sig, -log10(use.gene_df$FDR),decreasing = TRUE), ]

use.gene_df$gene   <- rownames(use.gene_df)

use.gene_df<- do.call(rbind, lapply(split(use.gene_df, use.gene_df$Sig), function(x) head(x, 100))) 

use.gene <- use.gene_df$gene

#mat <- TCGA_NSCLC_exp[use.gene,!grepl("11$", substr(colnames(TCGA_NSCLC_exp), 14, 15))] %>%as.matrix()
mat <- TCGA_NSCLC_exp[use.gene,] %>%as.matrix()
dim(mat)
cir1<-t(scale(t(mat)))
range(cir1)
mycol=colorRamp2(c(-2.5, 0.3, 3.1),c("blue","white","red"))#设置legend颜色，范围；可从https://www.58pic.com/peisebiao/网站进行配色

mycol1 = colorRamp2(c(-2, 0, 2), c("#003399", "white", "#cccc00"))

pdf('./Heatmap.pdf')
circos.heatmap(cir1,col=mycol1,dend.side = "inside",#聚类放在环形内测
               rownames.side = "outside")

dev.off()


circos.clear()#绘制完成后需要使用此函数完全清除布局













