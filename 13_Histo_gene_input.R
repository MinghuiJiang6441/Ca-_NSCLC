library(tidyr)
library(stringr)

Histo_gene <-read.csv('./data/Cytoscape_network1.csv')

head(Histo_gene)
Histo_gene$Samples <- str_sub(Histo_gene$Samples, 1, -3)

Histo_gene_long <- Histo_gene %>%
  pivot_longer(cols = -Samples, names_to = "Gene", values_to = "Expression")

head(Histo_gene_long)

top5_samples <- Histo_gene_long %>%
  group_by(Gene) %>%
  top_n(5, wt = Expression) %>%
  arrange(Gene, desc(Expression))

head(top5_samples)

write.csv(top5_samples, "./result/data/13_Histo_top5_samples.csv",quote = F, row.names = FALSE)

table(top5_samples$Samples)
