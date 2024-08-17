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
library(GSVA)
library(paletteer)
library(survival)
library(survminer)
library(cluster)
library(sva)
#

TCGA_NSCLC_exp<-readRDS("./data/TCGA_NSCLC_exp.rds")


Ca_Signatures <- unique(c(
  "MCU", "MICU1", "MICU2", "MICU3", "EMRE", "MCUb", 
  "MCUR1", "H-RAS12V", "FATE1", "STAT3", "AKT", 
  "PML", "p53", "TMX1", "PTEN","MCU", "EMRE", "MICU1", "MICU2", "MICU3", "MCUb", "MCUR1", "SLC25A23",
  "IP3R", "RyR", "SERCA", "VDAC", "mNCX", "NCLX", "FOXD1", "CREB",  "STIM1", "NPAS4", "CaMK",
  "miR-25", "miR-138", "miR-340", "EGR1", "SP1",  "PRMT1",
  "UCP2", "UCP3"
))

Ca_Signatures <- c("MICU1", "MICU2", "MICU3", "MCU", "MCUb", "EMRE", "MCUR1", "SLC25A23")



MCU_Geneset <- intersect(rownames(expr),Ca_Signatures)
MCU_Geneset
#[1] "MCUR1"    "MICU1"    "EGR1"     "SLC25A23" "PRMT1"    "TMX1"     "PML"      "FATE1"    "MICU3"    "MCU"      "MICU2"   
#[12] "STIM1"    "STAT3"    "PTEN"     "NPAS4"    "UCP3"     "UCP2"     "SP1"      "FOXD1"   
length(MCU_Geneset)
set.seed(0805)
MCU_Geneset_S <- sample(MCU_Geneset,size = 13)
saveRDS(MCU_Geneset,file = './result/data/MCU_Geneset.rds')
############################################################
#
#-------------------------处理临床数据----------------------
#
############################################################
process_clinical_data <- function(file_path, use_columns) {
  # 读取数据
  clinical_data <- read.delim(file_path, header = TRUE)
  
  # 提取有用的列
  clinical_data <- clinical_data[, use_columns]
  
  clinical_data$A0_Samples <- substr(x = clinical_data$A0_Samples,start = 1,stop = 12)
  # 删除重复的样本名
  clinical_data <- clinical_data[!duplicated(clinical_data$A0_Samples), ]
  
  # 设置行名为样本名
  rownames(clinical_data) <- clinical_data$A0_Samples
  
  # 返回处理后的数据框
  return(clinical_data)
}
file_path <- './data/LUAD_merge.cli.txt'

use_columns <- c("A0_Samples", "A1_Age", "race", "A7_OS.time", "A7_OS.event", "A6_gender", "ajcc_pathologic_stage", "A3_pathologic_N", "A4_pathologic_T", "A5_pathologic_stage")

LUAD_clin <- process_clinical_data(file_path, use_columns)


file_path <- './data/LUSC_merge.cli.txt'

LUSC_clin <- process_clinical_data(file_path, use_columns)

colnames(LUAD_clin) %in% colnames(LUSC_clin) %>% any()

TCGA_NSCLC_clin <- rbind(LUAD_clin,LUSC_clin)

head(TCGA_NSCLC_clin)



expr <-TCGA_NSCLC_exp[,!grepl("11$", substr(colnames(TCGA_NSCLC_exp), 14, 15))]
colnames(expr) <- substr(x = colnames(expr),start = 1,stop = 12)


any(colnames(expr) %in% rownames(TCGA_NSCLC_clin))

use.sample <- intersect(colnames(expr) , rownames(TCGA_NSCLC_clin))


TCGA_NSCLC_clin <- TCGA_NSCLC_clin[use.sample,]

expr<-expr[,use.sample]

colnames(expr)[1:4]

setdiff( rownames(TCGA_NSCLC_clin),colnames(expr) )

expr <- expr[,rownames(TCGA_NSCLC_clin)]

saveRDS(TCGA_NSCLC_clin,file = './data/02_TCGA_NSCLC_clin.rds')

saveRDS(expr,file = './data/02_TCGA_expr.rds')

TCGA_NSCLC_clin <- readRDS('./data/02_TCGA_NSCLC_clin.rds')

expr <- readRDS( './data/02_TCGA_expr.rds')

############################################################
#
#------------------------------ssgsea------------------------
#
############################################################
MCU_Geneset_S

gsva_mat <- as.data.frame(t(gsva(as.matrix(expr),gset.idx.list = list(MCU_Geneset),method = 'ssgsea'))) #计算得分

ssgsea_score <- cbind(gsva_mat,TCGA_NSCLC_clin)

names(ssgsea_score)[1]<- 'MCU_Score'

tail(ssgsea_score)

ssgsea_score$A7_OS.time <-as.numeric(ssgsea_score$A7_OS.time)
#ssgsea_score$A7_OS.event <- as.factor(ssgsea_score$A7_OS.event )
table(ssgsea_score$A7_OS.event)
ssgsea_score$A7_OS.event <- ifelse(ssgsea_score$A7_OS.event == "Dead", 1, 0) 

ssgsea_score <- ssgsea_score[-which(is.na(ssgsea_score$A7_OS.time )),]

ssgsea_score <- ssgsea_score [-which(ssgsea_score$A7_OS.time ==0) ,]

#ssgsea_score$MCU_Score_group <- ifelse( ssgsea_score$MCU_Score >= median(ssgsea_score$MCU_Score),"High","Low")%>%as.factor() #以中位数分组


# 将MCU_Score分为高低评分组
ssgsea_score$MCU_Score_Group <- ifelse(ssgsea_score$MCU_Score > median(ssgsea_score$MCU_Score), "High", "Low")

fit <- survfit(Surv(A7_OS.time , A7_OS.event) ~ MCU_Score_Group, data = ssgsea_score)

#saveRDS(ssgsea_score,'./result/data/ssgsea_score.rds')
ssgsea_score <- readRDS('./result/data/03_ssgsea_score.rds')
fit
pdf('./result/figure/03_Fig3A_Survplot.pdf',onefile = F)

ggsurvplot(
  fit, 
  data = ssgsea_score, 
  size = 1,                 # change line size
  palette = 
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

dev.off()

fit <- survfit(Surv(A7_OS.time , A7_OS.event) ~ MCU_Score_Group, data = ssgsea_score)

cutoff<-surv_cutpoint(ssgsea_score, #数据集
                      time="A7_OS.time",#“ ”里写数据集时间变量
                      event="A7_OS.event",##“ ”里数据集结局变量名称
                      variables='MCU_Score'
)



data_for_clustering <- ssgsea_score[, c("MCU_Score")]
# 进行K-means聚类（假设分为2个聚类）
set.seed(123) # 设置随机种子以保证结果可重复
kmeans_result <- kmeans(data_for_clustering, centers = 2)
ssgsea_score$Cluster <- as.factor(kmeans_result$cluster)
table(ssgsea_score$Cluster)
fit <- survfit(Surv(A7_OS.time, A7_OS.event) ~ Cluster, data = ssgsea_score)

ggsurvplot(fit, 
           data = ssgsea_score, 
           pval = TRUE, 
           conf.int = TRUE, 
           risk.table = TRUE, 
           ggtheme = theme_minimal(), 
           xlab = "时间 (天)", 
           ylab = "生存概率",
           title = "聚类结果的生存曲线")


















