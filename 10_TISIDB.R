TISIDB_geneset <- read.table('./data/TISIDB_set.txt' ,header=T, sep="\t", check.names=F)

Prognostic_gene

chemokine <- TISIDB_geneset[TISIDB_geneset[,2] == "chemokine", 1]
immunoinhibitor <- TISIDB_geneset[TISIDB_geneset[,2] == "Immunoinhibitor", 1]
immunostimulator <- TISIDB_geneset[TISIDB_geneset[,2] == "Immunostimulator", 1]
MHC <- TISIDB_geneset[TISIDB_geneset[,2] == "MHC", 1]
receptor <- TISIDB_geneset[TISIDB_geneset[,2] == "receptor", 1]
TIL <- TISIDB_geneset[TISIDB_geneset[,2] == "TIL", 1]

head(TISIDB_geneset)

library(linkET)
k <- intersect(immunostimulator,rownames(expr))

cor_res <- correlate(t(expr[Prognostic_gene,]), t(expr[k,]),method = "spearman")

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
  scale_fill_gradient2(low='#67B26F', high='#F2AA9D',mid = 'white',
                       limit=c(-1,1),name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_text(size=8,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())

plot1

ggsave(paste0('./result/figure/10_Fig10E_','immunostimulator','.pdf'),plot1,width = 5,height = 7)
















