expr <- readRDS('./data/02_TCGA_expr.rds')
library(estimate)
ls("package:estimate") 

in.file <- system.file("extdata", "sample_input.txt", package="estimate")
out.file <- tempfile(pattern="estimate", fileext=".gct")

#read.table(input.f, header = TRUE, row.names = 1, sep = "\t", quote = "", stringsAsFactors = FALSE)
Esti_data <- write.table(x = expr[,CoxData$Sample],file = './data/Esti_data.txt',sep = "\t", quote = FALSE,row.names = T)

#a <- read.table('./data/Esti_data.txt', header = TRUE, row.names = 1, sep = "\t", quote = "", stringsAsFactors = FALSE)

#rownames(a)[1:3]

filterCommonGenes(input.f='./data/Esti_data.txt',
                  output.f="./data/NSCLC_genes.gct",
                  id="GeneSymbol")

estimateScore(input.ds = "./data/NSCLC_genes.gct",
              output.ds="./data/NSCLC_estimate_score.gct", 
              platform="affymetrix")
# platform = c("affymetrix", "agilent", "illumina"))

plotPurity(scores="./data/NSCLC_estimate_score.gct", 
           platform="affymetrix",
           output.dir="./result/data/estimated_purity_plots")
#生成plot文件
scores=read.table("./data/NSCLC_estimate_score.gct",skip = 2,header = T)

rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])

head(scores,6)

scores<- as.data.frame(scores)
scores$SampleID <- rownames(scores)
scores$group <- CoxData$RiskGroup
save(scores,file = './result/data/estimate_score.rdata') 


scores$group<-factor(scores$group,levels=c('Low_Risk','High_Risk'))

vln_Plot <- list()
for(i in colnames(scores)[1:4]){
  p <- ggplot(data = scores,aes_string(x = 'group', 
                           y = i , 
                           fill = 'group'))+ 
    scale_fill_manual(values = c('#7EA6BF','#A64B5B')) +
    geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
                size = 0.8, color="black") +
    geom_boxplot(notch = TRUE, outlier.size = -1, 
                 color="black", lwd=0.8, alpha = 0.7) +
    geom_point(shape = 21, size=2, 
               position = position_jitterdodge(), 
               color="black", alpha = 0.6) +
    theme_bw() + 
    ylab(i) +
    xlab( 'Risk Group') +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))+
    stat_compare_means(method = "wilcox.test", label = "p.signif")  # 添加 Wilcoxon 检验结果
  vln_Plot[[i]] <- p
}

pdf('./result/figure/11_Fig11C_violinplot.pdf',width = 5,height = 6,onefile = F)
grid.arrange(grobs = vln_Plot, ncol = 2)
dev.off()









colnames(scores)
head(CoxData)

data1 <- CoxData[,c(Prognostic_gene,'RiskScore')]



cor_res <- correlate(data1,scores[1:4],method = "spearman")

# 先整理下数据
df_r <- cor_res$r %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene1") %>% 
  pivot_longer(-1,names_to = "gene2",values_to = "correlation")

df_p <- cor_res$p %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene1") %>% 
  pivot_longer(-1,names_to = "gene2",values_to = "pvalue")

df_cor <- df_r %>% 
  left_join(df_p) %>% 
  mutate(stars = cut(pvalue,breaks = c(-Inf,0.05,0.01,0.001,Inf),right = F,labels = c("***","**","*"," ")))
## Joining with `by = join_by(gene, cell_type)`

head(df_cor)

plot1<- ggplot(df_cor, aes(gene1,gene2))+
  geom_tile(aes(fill=correlation))+
  geom_text(aes(label=stars), color="black", size=4)+
  scale_fill_gradient2(low='#171559', high='#D94169',mid = 'white',
                       limit=c(-1,1),name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_text(size=8,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())

plot1

ggsave(filename = './result/figure/11_Fig11D_risk_TME_score.pdf',plot = plot1,width = 8,height = 4)

############################################################
#
#--------------TIDE score————————————————————————————————
#
############################################################

#TIDE数据标准化
TIDE_data <- t(apply(expr, 1, function(x)x-(mean(x)))) 

colnames(TIDE_data) <- seq_len(ncol(TIDE_data))


write.table(x = TIDE_data,file = './data/TIDE_data.txt',sep = "\t", quote = FALSE,row.names = T)






