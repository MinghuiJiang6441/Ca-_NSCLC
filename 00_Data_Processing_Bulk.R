library(BiocManager)
BiocManager::install("RTCGA")
BiocManager::install("RTCGA.mRNA")
library(RTCGA)
library(RTCGAToolbox)
library(RTCGA.mRNA)
library(GEOquery)
library(dplyr)
library(tidyr)
#RTCGA包用于获取TCGA数据
############################################################
#
#--------------TCGA数据下载—————————————————————————————————
#
############################################################

library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)

# 准备数据
TCGA-LUAD-seq <- GDCprepare(query)

query <- GDCquery(
  project = "TCGA-LUSC", # TCGA-LUSC 是肺鳞癌的项目名称
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)

getwd()
############################################################
#
#-------------------验证集：GSE19188———————————————————————
#
############################################################
# 加载必要的包
library(readr)



con <- gzfile("./data/GSE19188_series_matrix.txt.gz", "rt")
GSE19188 <- read.table(con,header = T,comment.char = "!",row.names = 1)

# 关闭连接
close(con)

GSE19188=as.matrix(GSE19188)


GSE_data_expr = as.data.frame(GSE19188) 

GSE_data_expr = data.frame(probe_id = rownames(GSE_data_expr),GSE_data_expr) 

GSE_gpl = read_csv(file = './data/GPL570.csv')

GSE_gpl = GSE_gpl[-which(grepl('///',GSE_gpl$`Gene Symbol`)),]

#names(GSE_gpl)[1] <- 'probe_id'

which(colnames(GSE_gpl) == 'Gene Symbol')

GSE_gpl = GSE_gpl[,c(1,11)]

head(GSE_gpl,3)

saveRDS(GSE_gpl,file = './data/GSE_gpl.rds')

probe_annotation = function(matrix,annotate,mathod = c('sum','mean','median','min','max')[2]){
  # matrix是一个表达矩阵，第一列为探针ID，其他列为每个探针ID对应样本的表达值
  # annotate是探针注释信息，包含两列吗，第一列为探针ID，第二列为探针ID的注释信息
  # mathod多个探针ID对应同一个symbol的处理方法，默认为均值
 
  names(matrix)[1] = 'probe_id'
  names(annotate) = c('probe_id','GeneSymbol')
  if(length(matrix$probe_id)==length(unique(matrix$probe_id))){
    if(length(annotate$probe_id)==length(unique(annotate$probe_id))){
      dat = merge(x = annotate,y = matrix) %>% dplyr::select(-probe_id)
      if(mathod == 'sum'){
        dat  =  aggregate(dat[,-1], by=list(symbol=dat[,1]),sum)
      }else if(mathod == 'mean'){
        dat  =  aggregate(dat[,-1], by=list(symbol=dat[,1]),mean)
      }else if(mathod == 'min'){
        dat  =  aggregate(dat[,-1], by=list(symbol=dat[,1]),min)
      }else if(mathod == 'median'){
        dat  =  aggregate(dat[,-1], by=list(symbol=dat[,1]),median)
      }else {
        dat  =  aggregate(dat[,-1], by=list(symbol=dat[,1]),max)
      }
      dat = as.data.frame(dat)
      return(dat)
    }else {
      print('输入的探针注释的probe ID有重复，请重新输入去重之后的探针注释文件')
    }
  }else {
    print('输入的探针表达矩阵中的probe ID有重复，请重新输入去重之后的探针表达矩阵')
  }
}

GSE_data_expr <- as.matrix(GSE_data_expr)
GSE_data_expr[, -1] <- apply(GSE_data_expr[, -1], 2, as.numeric)


symbol_exp = probe_annotation(matrix = GSE_data_expr,annotate = GSE_gpl,mathod = 'mean')

 rownames(symbol_exp) <-symbol_exp$symbol

GSE19188_exp <- symbol_exp[,-1]
GSE19188_exp[1:4,1:2]
saveRDS(GSE19188_exp,file = './data/GSE19188_exp.rds')


############################################################
#
#-------------------验证集：GSEGSE37745———————————————————————
#
############################################################

con <- gzfile("./data/GSE37745_series_matrix.txt.gz", "rt")
GSE37745 <- read.table(con,header = T,comment.char = "!",row.names = 1)

# 关闭连接
close(con)

GSE37745=as.matrix(GSE37745)

GSE_data_expr = as.data.frame(GSE37745) 
GSE_data_expr = data.frame(probe_id = rownames(GSE_data_expr),GSE_data_expr) 

symbol_exp = probe_annotation(matrix = GSE_data_expr,annotate = GSE_gpl,mathod = 'mean')

rownames(symbol_exp) <-symbol_exp$symbol

GSE37745_exp <- symbol_exp[,-1]

saveRDS(GSE37745_exp,file = './data/GSE37745_exp.rds')

library(readxl)
table_2 <- read.csv(file = "./data/clindata.csv",header = T)

head(table_2)
GSE19188_clin <- table_2

rownames(GSE19188_clin) <- GSE19188_clin$Accession

head(GSE19188_clin)
any(rownames(GSE19188_clin) %in% colnames(GSE19188_exp))

saveRDS(GSE19188_clin,'./data/GSE19188_clin.rds')



GSE37745_clin <- read.csv(file = "./data/GSE37745.ccsv",header = T)

head(GSE37745_clin)

rownames(GSE37745_clin) <- GSE37745_clin$Accession

any(rownames(GSE37745_clin) %in% colnames(GSE37745_exp))


saveRDS(GSE37745_clin,file = './data/GSE37745_clin.rds')



  