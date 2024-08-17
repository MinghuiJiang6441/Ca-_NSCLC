setwd("/data/home/jiangminghui/NSCLC")
library(glmnet)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(psych)
library(Hmisc)
library(GSVA) # BiocManager::install('GSVA')
library(GSEABase)
library(enrichplot)
library(limma)
library(ggplot2) 
library(stringr)
library(circlize)
library(ComplexHeatmap)
expr <- readRDS('./data/02_TCGA_expr.rds')

Prognostic_gene <- readRDS('./result/Prognostic_gene.rds')
Prognostic_gene

  cox_data <- readRDS('./result/data/CoxData.rds')
  ############################################################
  #
  #--------------------计算相关性矩阵————————————————————————————
  #
  ############################################################
head(cox_data,2)


Prognostic_gene

#saveRDS(Prognostic_gene,'./result/Prognostic_gene.rds')

#saveRDS(CoxData,file = './result/data/CoxData.rds')

gene_cor <- corr.test(t(df1[,1:2]), method = "pearson")
df1 <- expr[Prognostic_gene,] %>%as.matrix() 
df2 <- expr[sample(1:nrow(expr),50),]

correlation_results <- corr.test(t(df1), t(df2), method = "spearman")

result.r <- correlation_results$r  # 相关性矩阵
result.p <- correlation_results$p  # p值矩阵
saveRDS(cor.result,file ='./result/data/09_cor.result.rds')



row_means <- rowMeans( corrlist$r )
data_clean <- corrlist$r[, colSums(is.na(corrlist$r)) == 0]



############################################################
#
#---------------------------GSEA————————————————————————————
#
############################################################

#x <- readLines("./data/c2.cp.kegg.v7.4.symbols.gmt")
#res <- strsplit(x, "\t")
#names(res) <- vapply(res, function(y) y[1], character(1))
#genesets  <- lapply(res, "[", -c(1:2))
KEGGgmt <- read.gmt("./data/c2.cp.kegg.v7.4.symbols.gmt")

# 对基因进行排序
row_means <- sort(row_means, decreasing = TRUE)

egmt <- GSEA(row_means, TERM2GENE= KEGGgmt, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)

gseaplot2(egmt,geneSetID = rownames(gsea_results_df)[1],pvalue_table=T)
gseaplot2(egmt,geneSetID = rownames(gsea_results_df)[2],pvalue_table=T)
gseaplot2(egmt,geneSetID = rownames(gsea_results_df)[3],pvalue_table=T)


dotplot(egmt,showCategory=10)

class(genesets$KEGG_N_GLYCAN_BIOSYNTHESIS)

############################################################
#
#---------------------------Hallmark————————————————————————
#
############################################################

x <- readLines("./data/h.all.v2023.1.Hs.symbols.gmt")
res <- strsplit(x, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
genesets  <- lapply(res, "[", -c(1:2))


# 进行GSVA分析并保存结果
gsva_Hallmarker <- gsva(as.matrix(expr[, cox_data$Sample]), gset.idx.list = genesets, kcdf = "Gaussian", method = "gsva", parallel.sz = 5)
saveRDS(gsva_Hallmarker, file = "./result/data/gsva_Hallmarker.rds")

gsva_Hallmarker <-readRDS("./result/data/gsva_Hallmarker.rds")
# 构建设计矩阵和比较矩阵
group_list2 <- cox_data$RiskGroup
design <- model.matrix(~0 + factor(group_list2))
colnames(design) <- levels(factor(group_list2))
rownames(design) <- colnames(expr[, cox_data$Sample])
contrast.matrix <- makeContrasts(contrasts = "High_Risk-Low_Risk", levels = design)

# 线性拟合和贝叶斯检验
fit <- lmFit(gsva_Hallmarker, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff2 <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

# 处理结果数据
diff2$pathway <- rownames(diff2)
diff2$Program <- "-PosVSNeg"
result <- diff2
result$pathway=gsub("HALLMARK_","",result$pathway)
result$pathway=str_to_title(result$pathway)

k=0
gs=list()
data=data.frame()

k=k+1

df1=result
img<-df1[order(df1$t),]
row.names(img)=img$pathway

df2 = data.frame(x = factor(1:nrow(img),labels =row.names(img)), y = img$t,z=img$adj.P.Val)
str(df2)
df2 <- transform(df2, judge = ifelse(y>0, 'Yes', 'No'))
df2$judge=as.character(df2$judge)
df2$judge[df2$z>0.05] = "NS"
df2$judge=factor(df2$judge,levels=c('Yes','No',"NS"))
data=rbind(data,df2)

head(data,2)
df2_filtered <- df2

# 按 y 值排序并选取前五个 "Yes" 和 "No" 条目
df2_yes <- head(df2_filtered[df2_filtered$judge == "Yes", ], 5)
df2_no <- head(df2_filtered[df2_filtered$judge == "NS", ], 5)

# 合并筛选后的数据
df2_top10 <- rbind(df2_yes, df2_no)

pdf('./result/figure/09_Fig9B_Hallmark.pdf',width = 6,height = 4)
ggplot(data = df2_top10, mapping = aes(x = x, y = y, fill = judge)) + 
  geom_bar(stat = 'identity', position = 'identity') + 
  labs(title = "GSVA", x = "", y = "t value of GSVA score") +
  scale_fill_manual(values = c("#008020", "#08519C", "grey"), guide = FALSE) + 
  coord_flip() +
  theme_minimal(base_size = 15) +  # 使用简洁的主题
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  # 背景透明
    plot.background = element_rect(fill = "transparent", color = NA),   # 背景透明
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    axis.text = element_text(color = "black"),  # 坐标轴文字颜色
    axis.title = element_text(color = "black"),  # 坐标轴标题颜色
    plot.title = element_text(hjust = 0.5)  # 标题居中
  )

dev.off()

############################################################
#
#---------------------------Hallmark————————————————————————
#
############################################################

rownames(gsva_Hallmarker)=gsub("HALLMARK_","",rownames(gsva_Hallmarker))
rownames(gsva_Hallmarker)=str_to_title(rownames(gsva_Hallmarker))


top5_activated <- result %>% filter(t > 0) %>% arrange(desc(t)) %>% head(5)
top5_inhibited <- result %>% filter(t < 0) %>% arrange(t) %>% head(5)

top5_pathways <- rbind(top5_activated, top5_inhibited)

top5_pathway_names <- top5_pathways$pathway
# 重新排序样本，使High_Risk在前，Low_Risk在后
ordered_samples <- cox_data$Sample[order(cox_data$RiskGroup, decreasing = TRUE)]
gsva_Hallmarker_ordered <- gsva_Hallmarker[top5_pathway_names, ordered_samples]

# 创建列注释
column_ha <- HeatmapAnnotation(group = cox_data$RiskGroup[order(cox_data$RiskGroup, decreasing = TRUE)],
                               col = list(foo1 = c("High_Risk" = "#009491", "Low_Risk" = "#e9e29c")))

# 使用colorRamp2函数生成颜色映射函数
color_set <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
pdf('./result/figure/09_Fig9C_EnrichmentHeatmap.pdf',height = 6,width = 8)
# 绘制热图
Heatmap(gsva_Hallmarker_ordered, show_row_names = TRUE, 
        top_annotation = column_ha,
        show_column_names = FALSE, cluster_columns = FALSE, col = color_set)

dev.off()

############################################################
#
#---------------------------Hallmark————————————————————————
#
############################################################














