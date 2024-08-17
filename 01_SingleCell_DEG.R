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

GSE117570_SC <- readRDS(file = './data/GSE117570_SC.rds')

DimPlot(GSE117570_SC, reduction = "umap", group.by = "seurat_clusters",cols = random_colors)

Idents(GSE117570_SC). <- 'seurat_clusters'
GSE117570_SC <- JoinLayers(GSE117570_SC)


GSE117570_SCe_marke <- FindAllMarkers(GSE117570_SC,only.pos = T,min.pct = 0.25, 
                                      logfc.threshold = 1,
                                      min.diff.pct = 0.25)

GSE117570_SC <- FindClusters(GSE117570_SC, resolution = 0.8)

saveRDS(GSE117570_SC, file = "./data/GSE117570_SC.rds")
############################################################
#
#--------------细胞注释————————————————————————————————
#
############################################################

ref <- celldex::HumanPrimaryCellAtlasData()
seurat_matrix <- GetAssayData(GSE117570_SC, layer = "data")

singleR_results <- SingleR(test = seurat_matrix, ref = ref, labels = ref$label.main)
GSE117570_SC <- AddMetaData(GSE117570_SC, metadata = singleR_results$labels, col.name = "SingleR_Labels")
GSE117570_SC

#设置颜色
d_palettes<- palettes_d_names

mycol<-paletteer_d( "ggsci::default_igv",n=51)
set.seed(123) # 为了可重复性设置随机种子
random_colors <- sample(mycol, 12)


DimPlot(GSE117570_SC, reduction = "umap",label = TRUE,
        group.by = c("SingleR_Labels",'seurat_clusters'),
        cols = sample(mycol))
class(ref)



Neutrophils <- c('CXCR2','FCGR3B')
Macrophages <- c('CD68', 'CD163', 'CD14')
Tcells <- c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','TNFRSF18','TNFRSF25')
Endothelial  <-c('PECAM1', 'VWF')
NK <-c('KLRB1','NCR1')
Fibroblast <-c('FGF7','MME', 'ACTA2','VCAN')
Bcells <-c('CD19', 'CD79A', 'MS4A1')

markergeneList <- list(Neutrophils=Neutrophils,Macrophages=Macrophages,Tcells=Tcells,
                      Endothelial=Endothelial,NK=NK,Fibroblast=Fibroblast, 
                    Bcells=Bcells)


Epithelial <- c(10,6,8,13,12,4)
Macrophages <- c(7)
Monocyte<- c(0,11,3,14)
NK <- c(1)
Tcells <- c(2,5)
Bcells <- c(9)
setdiff(0:14,c(Epithelial,Macrophages,Monocyte,NK,Tcells,Bcells))

current.cluster.ids <- c(Epithelial,Macrophages,Monocyte,NK,Tcells,Bcells)

new.cluster.ids <- c(rep("Epithelial",length(Epithelial)),
                     rep("Macrophages",length(Macrophages)),
                     rep("Monocyte",length(Monocyte)),
                     rep("NK",length(NK)),
                     rep("Tcells",length(Tcells)),
                     rep("Bcells",length(Bcells)))

GSE117570_SC@meta.data$Subcluster <- plyr::mapvalues(x = GSE117570_SC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
head(GSE117570_SC@meta.data)

pdf(  './result/figure/00_Fig1E_UMAP.pdf',width = 7,height = 6 )
DimPlot(GSE117570_SC, reduction = "umap", group.by = "seurat_clusters",cols = mycol)
dev.off()

pdf('./result/01_Fig1G_UMAP_annotation.pdf',width = 7,height = 6)
DimPlot(GSE117570_SC, reduction = "umap", group.by = "Subcluster",cols = random_colors)
dev.off()

############################################################
#
#--------------------差异基因鉴定————————————————————————————
#
############################################################


# 筛选每个细胞簇的标记基因
Idents(GSE117570_SC) <-"Subcluster" 
markers <- FindAllMarkers(GSE117570_SC, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 1, 
                          test.use = "wilcox")

# 过滤标记基因
markers <- markers[markers$p_val_adj < 0.05, ]

head(markers)
# 选取每组前5个基因
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

  # 画气泡图
pdf('./result/01_Fig1F_DotPlot.pdf',width = 7,height = 6)
DotPlot(GSE117570_SC, features = unique(top_markers$gene) )+ RotatedAxis()+ 
  scale_color_gradient(low = "white", high = "red")
dev.off()


# 根据orig.ident添加标签
GSE117570_SC$group <- ifelse(grepl("Tumor", GSE117570_SC$orig.ident), "NSCLC", "Control")

Idents(GSE117570_SC) <- 'group'

diff_cells <- FindMarkers(GSE117570_SC, 
                          ident.1 = "NSCLC", 
                          ident.2 = "Control",
                          test.use = "wilcox", 
                          min.pct = 0.25, logfc.threshold = 0)


Idents(GSE117570_SC) <- 'group'
# 提取NSCLC和Control的基因名

nsclc_genes <- rownames(diff_cells[diff_cells$avg_log2FC > 0, ])
control_genes <- rownames(diff_cells[diff_cells$avg_log2FC < 0, ])

# 提取表达矩阵
expr_matrix <- as.matrix(GetAssayData(GSE117570_SC, layer = "data"))

# 定义基因集
gene_sets <- list(
  NSCLC = nsclc_genes,
  Control = control_genes
)

# 运行ssGSEA
ssgsea_scores <- gsva(expr_matrix, gene_sets, method = "ssgsea")

# 将ssGSEA得分添加到Seurat对象中
GSE117570_SC[["ssgsea"]] <- CreateAssayObject(data = ssgsea_scores)

saveRDS(GSE117570_SC, file = "./data/GSE117570_SC.rds")


Idents(GSE117570_SC)<-'seurat_clusters'
markers <- FindMarkers(object =GSE117570_SC,group.by = 'group',  
                       logfc.threshold = 1,ident.1 = 'NSCLC',ident.2 = 'Control',
                       only.pos = FALSE,
                       test.use = "wilcox")
dim(markers) #1730    5


# 筛选符合阈值的基因
significant_markers <- markers[abs(markers$avg_log2FC) >= 1 & markers$p_val_adj < 0.05, ]

dim(significant_markers) # 1511


save(significant_markers,file = './result/data/01_significant_markers.RData')






# 识别细胞簇
Idents(GSE117570_SC) <- "Subcluster"

DefaultAssay(GSE117570_SC) <-'ssgsea'
# 提取NSCLC和对照组的ssGSEA得分
nsclc_scores <- FetchData(GSE117570_SC, vars = "ssgsea_NSCLC")
control_scores <- FetchData(GSE117570_SC, vars = "ssgsea_Control")

# 创建一个数据框来存储这些得分
scores_df <- data.frame(
  cell = colnames(GSE117570_SC),
  nsclc_scores = nsclc_scores,
  control_scores = control_scores,
  subcluster = GSE117570_SC@meta.data$Subcluster
)
scores_long <- scores_df %>%
  pivot_longer(cols = c(ssgsea_NSCLC, ssgsea_Control), 
               names_to = "group", 
               values_to = "score")

saveRDS(scores_long ,'./data/01_NSCLC_Control_ssGSEA.rds')


p <- ggplot(scores_long, aes(x = subcluster, y = score, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  theme_minimal() +
  labs(title = "NSCLC vs Control ssGSEA Scores by Subcluster",
       x = "Subcluster",
       y = "ssGSEA Score") +
  stat_compare_means(aes(group = group), method = "t.test", label = "p.signif")+
  scale_fill_manual(values = c("ssgsea_NSCLC" = "#FF9999", "ssgsea_Control" = "#99CCFF"))

pdf('./result/01_Fig2A_NSCLC_Control_ssGSEA_boxplot.pdf',width = 16,height = 9)
p
dev.off()


Idents(GSE117570_SC) <-'seurat_clusters'
DefaultAssay(GSE117570_SC) <-'RNA'

#markers <- FindAllMarkers(GSE117570_SC, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)



# Filter markers based on the given thresholds
significant_markers <- markers %>%
  filter(abs(avg_log2FC) >= 1 & p_val_adj < 0.05)
significant_markers$group <- ifelse(significant_markers$avg_log2FC >0,'pos','nag')

# Select the top 5 genes for each cluster based on the absolute value of avg_log2FC
top_genes <- significant_markers %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 10) %>%
  ungroup()

saveRDS(c(significant_markers,top_genes),file = './result/data/01_Fig2B_ManhattanPlot.rds')


# Create a scatter plot with jitter
p <- ggplot(significant_markers, aes(x = cluster, y = avg_log2FC, color = avg_log2FC > 0)) +
  geom_jitter(size = 3, alpha = 0.7, width = 0.2) +
  scale_color_manual(values = c("TRUE" = "#B70000", "FALSE" = "#0092B2")) +
  theme_minimal() +
  labs(title = "Differential Gene Expression by Cluster",
       x = "Cluster",
       y = "Average Log2 Fold Change",
       color = "Log2 Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_text(data = top_genes, aes(label = gene), vjust = 1.5, size = 3, check_overlap = TRUE, color = "black")+
  geom_text(aes(label = cluster), y = 0, vjust = 0.5, size = 4, color = "black")


pdf('./result/figure/01_Fig2B_ManhattanPlot.pdf',width = 16,height = 9)
p
dev.off()

head(significant_markers)
binary_matrix <- significant_markers %>%
  select(gene, cluster) %>%
  mutate(value = 0) %>%
  pivot_wider(names_from = cluster, values_from = value, values_fill = list(value = 1)) %>%
  rename_with(~paste0("cluster", .), -gene)

# 将其他列添加回去
final_data <- significant_markers %>%
  select(-cluster) %>%
  distinct() %>%
  left_join(binary_matrix, by = "gene")
head(final_data)

saveRDS(c(significant_markers,final_data,binary_matrix),file ='./result/data/01_Fig2C_Upset_Plot.rds' )
# 绘制UpSet图
pdf('./result/figure/01_Fig2C_Upset_Plot.pdf',width = 10,height = 8)
upset(final_data)

dev.off()
head(binary_df)

DEGs1<- unique(significant_markers$gene)
saveRDS(DEGs1,file = './result/data/01_DEGs1.rds')





significant_markers <-readRDS('./result/data/significant_markers1.rds')

significant_markers <- significant_markers %>%
  filter(abs(avg_log2FC) >= 1 & p_val_adj < 0.05)
head(significant_markers)
# 添加一个唯一标识符列
significant_markers <- significant_markers %>%
  mutate(id = row_number())

# 将significant_markers转化为列为不同cluster，基因竖着排列
significant_markers_wide <- significant_markers %>%
  dplyr::select(id, cluster, gene) %>%
  tidyr::pivot_wider(names_from = cluster, values_from = gene)

# 将结果写入txt文件

na.omit(significant_markers_wide[2])

l <-split(significant_markers$gene,significant_markers$cluster)
a <-do.call(cbind,l)
write.table(a, file = "./result/data/01_significant_markers_by_cluster.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

dim(a)
# 打开一个txt文件进行写入
fileConn <- file("./result/data/01_significant_markers_by_cluster.txt")

pdf('./result/figure/01_Fig2C_Upset_Plot.pdf',onefile = F)
upset(fromList(l),order.by = "freq")
dev.off()
