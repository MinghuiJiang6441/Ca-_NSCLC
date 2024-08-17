setwd("/data/home/jiangminghui/NSCLC")
library(preprocessCore)

library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
#BiocManager::install("RcolorBrewer")
#install.packages("RColorBrewer")
library(RColorBrewer)
library(tidyHeatmap)
library(tidyverse)
library(RColorBrewer)
library(tidyr)
library(tibble)
library(ggsci)
library(IOBR)

source('./data/CIBERSORT.R')

CoxData <-readRDS('./result/data/CoxData.rds')
#LM22.file <- "./data/LM22.txt"

LM22.file<- read.table(file = "./data/LM22.txt", header=T, sep="\t", row.names=1,check.names=F)


GEO_cibersort.results <- CIBERSORT(LM22.file ,expr, perm = 50, QN = T)
#saveRDS(GEO_cibersort.results,'./result/data/GEO_cibersort.results.rds')
GEO_cibersort.results <-readRDS('./result/data/GEO_cibersort.results.rds')
GEO_cibersort.results <- GEO_cibersort.results[rownames(CoxData),]


cibersort_data <- as.data.frame(GEO_cibersort.results[,1:22])
cibersort_data<-rownames_to_column(cibersort_data,var="Sample")

Group<- CoxData$RiskGroup

cibersort<-cbind(cibersort_data,Group)
cibersort<- melt(cibersort,id.vars=c("Sample","Group"))

colnames(cibersort)<-c("Sample","Group","celltype","composition") 
#连接分组信息和cibersort结果文件

#绘制差异分析箱线图
boxplot_cibersort<- ggplot(cibersort, aes(x = celltype, y = composition))+ 
  labs(y="Cell composition",x= "")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5)+ 
  scale_fill_npg()+
  #修改主题
  theme_bw() + 
  theme(axis.title = element_text(size = 12,color ="black"), 
        axis.text = element_text(size= 12,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)
  ) +
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)#隐藏不显著的

ggsave(file="./result/figure/10_Fig10C_immunecellshigh-lowgroups.pdf",boxplot_cibersort,height=10,width=15)

ggplot(cibersort, aes(x = celltype, y = composition, fill = Sample)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Group", y = "Relative Abundance", fill = "Immune Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

res_cibersort <- cell_bar_plot(GEO_cibersort.results)


p1 <- cibersort %>% 
  ggplot(aes(Sample,composition))+
  geom_bar(stat = "identity",position = "stack",aes(fill=celltype))+
  labs(x=NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = palette4,name=NULL)+ # iobr还给大家准备了几个色盘，贴心！
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
  )
p1
ggsave(file="./result/figure/10_Fig10A_Cibersort_Abundance.pdf",p1,height=10,width=15)



# 有顺序的箱线图
library(forcats)


p2 <- ggplot(cibersort,aes(fct_reorder(celltype, composition),composition,fill = celltype)) + 
  geom_boxplot() + 
  #geom_jitter(width = 0.2,aes(color=cell_type))+
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = palette4)
p2
ggsave(file="./result/figure/10_Sup_Fig_Cibersort_boxplot.pdf",p2,height=10,width=15)


library(ggpubr)
library(stringr)

# 分组
head(cibersort,2)

p3 <- ggplot(cibersort,aes(fct_reorder(celltype,composition),composition,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = Group,label = ..p.signif..),
                     method = "kruskal.test",label.y = 0.4)
p3
ggsave(file="./result/figure/10_Fig10C_immunecellshigh-lowgroups2.pdf",p3,height=10,width=15)


immune_data <- GEO_cibersort.results[, 1:22]

library(psych)

dim(immune_data)
# 计算Spearman相关性
cor_matrix <- corr.test(immune_data, method = "spearman")
# 加载包
library(corrplot)

# 把数据由数据框转换为矩阵格式
corr <- as.matrix(cor_matrix$r)

pvalue <- as.matrix(cor_matrix$p)

#设置颜色
addcol <- colorRampPalette(c("blue", "white", "red")) 

pdf('./result/figure/10_Fig10B_immune_cells_cor.pdf',width = 8,height = 8,onefile = F)
corrplot(corr, # 相关性矩阵 
         method = "color", # 表示用颜色表示相关性大小
         col = addcol(100), 
         tl.col = "black", # 文本标签的颜色
         tl.cex = 0.8, # 文本标签的字符大小
         tl.srt = 90, #  文本标签的旋转角度
         tl.pos = "td", # 文本标签位置，td表示顶部和对角线 
         p.mat = pvalue, #  P 值矩阵
         diag = T, # 是否显示对角线上的相关性值
         type = 'upper', # 只绘制上三角部分
         sig.level = c(0.05), # 设置显著性水平阈值，可设置多个
         pch.cex = 1,  # 显著性标记字符大小
         pch.col = 'grey20',  # 显著性标记字符颜色
         insig = 'label_sig',
         order = 'AOE', #设置一种排序方式
)
dev.off()


cor_matrix1 <- corr.test(x = immune_data,y = CoxData$RiskScore, method = "spearman")

df_r <- cor_matrix1$r %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "cell_type") %>% 
  pivot_longer(-1,names_to = "risk",values_to = "correlation")

df_p <- cor_matrix1$p %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "cell_type") %>% 
  pivot_longer(-1,names_to = "risk",values_to = "pvalue")

df_cor <- df_r %>% 
  left_join(df_p) %>% 
  mutate(stars = cut(pvalue,breaks = c(-Inf,0.05,0.01,0.001,Inf),right = F,labels = c("***","**","*"," ")))
## Joining with `by = join_by(gene, cell_type)`

head(df_cor)
# 以EGFR为例
df_egfr <- df_cor %>% 
  filter(cell_type1=="B cells naive")

text_color <- c(rep("red",4),rep("#FDD819",1),rep("#67B26F",1),rep("#C4E0E5",12),rep("#67B26F",2),rep("#FDD819",3),rep("red",5))

pdf(file = './result/figure/10_Fig10D_segment_plot.pdf',height = 5,width = 6)
p <-ggplot(df_cor, aes(correlation, fct_reorder(cell_type,correlation)))+
  geom_segment(aes(xend = 0,yend = cell_type),color="grey70",size=1.2)+ 
  geom_point(aes(size = abs(correlation), color=factor(stars)))+ 
  scale_color_manual(values = c("red","#FDD819","#67B26F","#C4E0E5"),name="P-value")+
  scale_size_continuous(range = c(3,8),name="Correlation")+ 
  labs(x=NULL,y=NULL)+
  theme_bw()+
  theme(axis.text.x = element_text(color="black",size = 14),
        axis.text.y = element_text(color = text_color,size = 14)#给标签上色没想到特别好的方法
  )
dev.off()
ggsave('./result/figure/10_Fig10D_segment_plot.pdf',p ,height = 6,width = 7)







