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

library(igraph)

library(RColorBrewer)

library(paletteer)
############################################################
#
#--------------GSE117570单细胞数据初始化处理————————————————

#--------------读数原始数据转换成Seurat对象------------------
#
############################################################

#解压zip格式文件
utils::unzip('./data/GSE117570_RAW.zip', exdir = "./data")#unzip

#将文件储存到list中
data_dir <- "./data/GSE117570_RAW"
files <- list.files(data_dir, full.names = TRUE, pattern = ".txt")

# 创建一个空列表来存储Seurat对象
scRNAlist <- list()

# 遍历每个文件，读取数据并创建Seurat对象
for (file in files) {
  # 读取txt文件数据
  counts <- read.table(file, header =TRUE, check.names = TRUE) %>% as.matrix # 假设第一列是基因名
  counts <- na.omit(counts)
  # 从文件名中提取样本名，假设文件名格式为"GSMxxxxxx_XX_SampleType_processed_data.txt"
  sample_name <- sub("_processed_data.txt", "", basename(file))
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = counts,
                                   assay = 'RNA',
                                   project = "GSE117570",
                                   min.cells = 3,
                                   min.features = 200)
 
  
  # 将样本名添加到元数据中
  seurat_obj@meta.data$orig.ident <- sample_name

  # 将Seurat对象添加到列表中
  scRNAlist[[sample_name]] <- seurat_obj
}


names(scRNAlist)

GSE117570_SC <-merge(x =scRNAlist[[1]], y = scRNAlist[-1] ) #合并样本

GSE117570_SC

# 查看合并和过滤后的Seurat对象的维度
print(dim(GSE117570_SC@assays$RNA@data))

# 打印合并后每个样本的原始标识符数量
print(table(GSE117570_SC@meta.data$orig.ident))


GSE117570_SC

############################################################
#
#------------------检测数据并初步过滤--------————————————————
#
############################################################


# 计算线粒体基因比例
GSE117570_SC[["percent.mt"]] <- PercentageFeatureSet(GSE117570_SC, pattern = "^MT-")

# 查看前6个细胞的QC指标
head(GSE117570_SC@meta.data, 6)

VlnPlot(GSE117570_SC,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"), #纵坐标
        ncol = 3
)

plot1 <- FeatureScatter(GSE117570_SC,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(GSE117570_SC,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1 + plot2

GSE117570_SC <- subset(GSE117570_SC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

GSE117570_SC <-SCTransform(object = GSE117570_SC,assay = "RNA",vars.to.regress = c('nFeature_RNA','nCount_RNA'),verbose = TRUE)
DefaultAssay(GSE117570_SC)
############################################################
#
#-------------------------数据标准化--------————————————————
#
############################################################
#标准化数据：
GSE117570_SC <- NormalizeData(GSE117570_SC, verbose = FALSE,scale.factor = 10000)
#GSE117570_SC <- RunALRA(GSE117570_SC)

GSE117570_SC[['RNA']]@data
############################################################
#
#------------------------筛选高变量基因--------———————————
#
############################################################

#筛选高变量基因"
  GSE117570_SC <- FindVariableFeatures(GSE117570_SC, nfeatures = 2000,
                                       assay = 'RNA',
                                       selection.method = "vst") 

top10 <- head(VariableFeatures(GSE117570_SC), 10)# 筛选前10个高变基因
top10
DefaultAssay(GSE117570_SC)  <- 'RNA'
plot1 <- VariableFeaturePlot(GSE117570_SC) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge = 0,ynudge = 0) # 额外标记前10个HVG。
plot1 + plot2 #图层叠加

pdf('./result/00_Fig1B_VariableFeaturePlot.pdf',width = 12,height = 6)
plot1 + plot2
dev.off()



  #二次标准化"

all.genes <- rownames(GSE117570_SC)

GSE117570_SC <- ScaleData(GSE117570_SC ,features = all.genes) 

GSE117570_SC <- ScaleData(GSE117570_SC) # 默认对HVG进行二次标准化

GSE117570_SC <- ScaleData(GSE117570_SC, vars.to.regress = "percent.mt")  
 


 #降维:
  GSE117570_SC <- RunPCA(GSE117570_SC,npcs = 10,features = VariableFeatures(object = GSE117570_SC)) 
  DimHeatmap(GSE117570_SC, dims = 1:10, cells = 500, balanced = TRUE)
  
  VizDimLoadings(GSE117570_SC, dims = 1:10, reduction = "pca")
  
  DimPlot(GSE117570_SC, reduction = "pca")
  
  GSE117570_SC <- JackStraw(GSE117570_SC, num.replicate = 100)
  GSE117570_SC <- ScoreJackStraw(GSE117570_SC, dims = 1:5)
  JackStrawPlot(GSE117570_SC,xmax = 0.01, ymax = 0.1,dims = 1:5) # 将重要的PC的P值进行可视化，通常情况下此步骤可以跳跃，选择不同数量的PC结果会有很大的不同，且建议使用更多的PC进行下游分析，如果选择过少，则会丢失较多的特征值，对分析产生不利影响。通常默认选择20个，可选择10~50个，结果通常产生不了太大变化。
  ElbowPlot(GSE117570_SC)
  
  pdf('./result/00_Fig1C_JackStrawPlot.pdf',width = 5,height = 5)
  JackStrawPlot(GSE117570_SC, dims = 1:5)
  dev.off()
  
  pdf('./result/00_Fig1D_ElbowPlot.pdf',width = 5,height = 5)
  ElbowPlot(GSE117570_SC)
  dev.off()
  

 #聚类
  
GSE117570_SC <- FindNeighbors(GSE117570_SC, dims = 1:10)

GSE117570_SC <- FindClusters(GSE117570_SC, resolution = 0.5)

GSE117570_SC <-  RunUMAP(GSE117570_SC,dims = 1:10) %>% RunTSNE(dims = 1:10)

GSE117570_SC <-RunHarmony(GSE117570_SC,reduction = 'pca',lambda = 1,group.by.vars = 'orig.ident')
#可视化：
#使用DimPlot函数可视化聚类结果：
DimPlot(GSE117570_SC, reduction = "umap", group.by = "seurat_clusters")


d_palettes<- palettes_d_names

mycol<-paletteer_d( "ggsci::default_igv",n=51)
set.seed(123) # 为了可重复性设置随机种子
random_colors <- sample(mycol, 12)

pdf('./result/00_Fig1E_UMAP.pdf',width = 8,height = 6)
DimPlot(GSE117570_SC, reduction = "umap", group.by = "seurat_clusters",cols = random_colors)
dev.off()


saveRDS(GSE117570_SC, file = "./data/GSE117570_SC.rds")



