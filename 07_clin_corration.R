library(dplyr)
library(pheatmap)
library(tibble)
head(TCGA_NSCLC_clin,2)
head(risk_data,2)

colnames(TCGA_NSCLC_clin)

use.sample <- rownames(risk_data)
ClinicalFeatures <- TCGA_NSCLC_clin[use.sample,]

ClinicalFeatures <- cbind(ClinicalFeatures,risk_data$RiskScore)

ClinicalFeatures <-ClinicalFeatures[,-11]

names(ClinicalFeatures)[11] <- 'RiskScore'

head(ClinicalFeatures,2)

names(ClinicalFeatures)[c(6,9,8,7)] <- c('gender','T.stage','N.stage','M.stage')

ClinicalFeatures$Age <- ifelse(ClinicalFeatures$A1_Age >= 60,'>=60','<60')
head(ClinicalFeatures,2)

table(ClinicalFeatures$race)

age_plot <-ggplot(ClinicalFeatures, aes(x = Age, y = RiskScore, fill = Age)) +
  geom_boxplot(aes(color = Age), fill = "white", size = 1.2) +  # 设置箱子边框加粗
  geom_jitter(aes(color = Age), size = 2) +  # 设置点的大小加大
  stat_compare_means(method = "anova", label.y = 6) +
  scale_color_manual(values = c("yellow", "#56B0F5")) +
  theme_bw() +
  ggtitle("Age")  # 添加主题
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )+
    geom_signif(comparisons = list(c("<60", ">=60")), 
                map_signif_level = TRUE)
  

  
 gender_plot <-  ggplot(ClinicalFeatures, aes(x = gender, y = RiskScore, fill = gender)) +
   geom_boxplot(aes(color = gender), fill = "white", size = 1.2) +  # 设置箱子边框加粗
   geom_jitter(aes(color = gender), size = 2) +  # 设置点的大小加大
   stat_compare_means(method = "anova", label.y = 6) +
   scale_color_manual(values = c("yellow", "#56B0F5")) +
   theme_bw() +
   ggtitle("gender") +  # 添加主题
   theme(
     plot.title = element_text(hjust = 0.5),
     legend.position = "none"
   ) +
   geom_signif(comparisons = list(c("male", "female")), 
               map_signif_level = TRUE)
  
 colnames(ClinicalFeatures) 
 filtered_data <- ClinicalFeatures %>%
   filter(race %in% c("black or african american", "white"))
 
 race_plot <-  ggplot(filtered_data, aes(x = race, y = RiskScore, fill = race)) +
   geom_boxplot(aes(color = race), fill = "white", size = 1.2) +  # 设置箱子边框加粗
   geom_jitter(aes(color = race), size = 2) +  # 设置点的大小加大
   stat_compare_means(method = "anova", label.y = 6) +
   scale_color_manual(values = c("yellow", "#56B0F5")) +
   theme_bw() +
   ggtitle("race") +  # 添加主题
   theme(
     plot.title = element_text(hjust = 0.5),
     legend.position = "none"
   ) +
   geom_signif(comparisons = list(c("black or african american", "white")), 
               map_signif_level = TRUE)
 
 
 
 race_plot <-  ggplot(filtered_data, aes(x = race, y = RiskScore, fill = race)) +
   geom_boxplot(aes(color = race), fill = "white", size = 1.2) +  # 设置箱子边框加粗
   geom_jitter(aes(color = race), size = 2) +  # 设置点的大小加大
   stat_compare_means(method = "anova", label.y = 6) +
   scale_color_manual(values = c("yellow", "#56B0F5")) +
   theme_bw() +
   ggtitle("race") +  # 添加主题
   theme(
     plot.title = element_text(hjust = 0.5),
     legend.position = "none"
   ) +
   geom_signif(comparisons = list(c("black or african american", "white")), 
               map_signif_level = TRUE)
 
 
 
 ClinicalFeatures <- ClinicalFeatures %>%
   mutate(T.stage = case_when(
     T.stage %in% c("T1", "T1a", "T1b", "T2", "T2a", "T2b") ~ "T1/2",
     T.stage %in% c("T3", "T4") ~ "T3/4",
     TRUE ~ T.stage
   ))
 
 
 table(ClinicalFeatures$T.stage)
 
 filtered_data <- ClinicalFeatures %>%
   filter(T.stage %in% c("T1/2", "T3/4"))

 T.stage_plot <- ggplot(filtered_data, aes(x = T.stage, y = RiskScore, fill = T.stage)) +
   geom_boxplot(aes(color = T.stage), fill = "white", size = 1.2) +  # 设置箱子边框加粗
   geom_jitter(aes(color = T.stage), size = 2) +  # 设置点的大小加大
   stat_compare_means(method = "anova", label.y = 6) +
   scale_color_manual(values = c("yellow", "#56B0F5")) +
   theme_bw() +
   ggtitle("T.stage") +  # 添加主题
   theme(
     plot.title = element_text(hjust = 0.5),
     legend.position = "none"
   ) +
   geom_signif(comparisons = list(c("T1/2", "T3/4")), 
               map_signif_level = TRUE) 
 
 table(ClinicalFeatures$N.stage)
 
 
 filtered_data <- ClinicalFeatures %>%
   mutate(N.stage = case_when(
     N.stage == "N0" ~ "N0",
     N.stage %in% c("N1", "N2") ~ "N1/2",
     TRUE ~ NA_character_
   )) %>%
   filter(!is.na(N.stage)) 
 
 N.stage_plot <-ggplot(filtered_data, aes(x = N.stage, y = RiskScore, fill = N.stage)) +
   geom_boxplot(aes(color = N.stage), fill = "white", size = 1.2) +  # 设置箱子边框加粗
   geom_jitter(aes(color = N.stage), size = 2) +  # 设置点的大小加大
   stat_compare_means(method = "anova", label.y = 6) +
   scale_color_manual(values = c("yellow", "#56B0F5")) +
   theme_bw() +
   ggtitle("N.stage") +  # 添加主题
   theme(
     plot.title = element_text(hjust = 0.5),
     legend.position = "none"
   ) +
   geom_signif(comparisons = list(c("N0", "N1/2")), 
               map_signif_level = TRUE) 
  


combined_plot <- grid.arrange(age_plot, gender_plot, N.stage_plot, race_plot, T.stage_plot, ncol = 2)

ggsave('./result/figure/07_Fig7B_combined_plot.pdf', combined_plot, width = 8, height = 12)

############################################################
#
#--------------------相关性热图———————————————————————————
#
###########################################################
library(Hmisc)
library(corrplot)#先加载包

cor_data <- risk_data[,c('RiskScore',Prognostic_gene)]
head(cor_data)

Hmisc::rcorr(as.matrix(cor_data), type = "pearson") -> corrlist

#options(repr.plot.width = 11, repr.plot.height=11)
pdf('./result/figure/07_Fig7C_Correlation_heatmap.pdf',width = 5,height = 5,onefile = FALSE)
corrplot(as.matrix(corrlist$r),
         col = rev(COL2('RdBu', 200)),
         type = "upper",
         p.mat = corrlist$p, 
         order = "hclust",
         sig.level = 0.01, 
         insig = "blank",
         tl.col = "black", 
         tl.srt = 60)

dev.off()


############################################################
#
#--------------------表达分布热图———————————————————————————
#
###########################################################

ClinicalFeatures <- TCGA_NSCLC_clin[use.sample,]

ClinicalFeatures <- cbind(ClinicalFeatures,risk_data[,-c(1,2)])

ClinicalFeatures$Age <- ifelse(ClinicalFeatures$A1_Age >= 60,'>=60','<60')

table(ClinicalFeatures$A3_pathologic_N)

colnames(ClinicalFeatures)

ClinicalFeatures <- ClinicalFeatures %>%
  mutate(A4_pathologic_T = case_when(
    A4_pathologic_T %in% c("T1", "T1a", "T1b") ~ "T1",
    A4_pathologic_T %in% c("T2", "T2a", "T2b") ~ "T2",
    A4_pathologic_T %in% c("T3") ~ "T3",
    A4_pathologic_T %in% c( "T4") ~ "T4",
    
    TRUE ~ A4_pathologic_T
  ))


ClinicalFeatures <- ClinicalFeatures %>%
  mutate(A5_pathologic_stage = case_when(
    A5_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "Stage 1",
    A5_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage 2",
    A5_pathologic_stage %in% c("Stage III",'Stage IIIA','Stage IIIB') ~ "Stage 3",
    A5_pathologic_stage %in% c( "Stage IV") ~ "Stage 4",
    
    TRUE ~ A5_pathologic_stage
  ))

ClinicalFeatures <- ClinicalFeatures %>%
  rename(
    Stage = A5_pathologic_stage,
    N.stage = A3_pathologic_N,
    T.stage = A4_pathologic_T,
    Sample = A0_Samples,
    OS.time = A7_OS.time,
    OS = A7_OS.event,
    gender = A6_gender
  )

ClinicalFeatures <- ClinicalFeatures %>%
  mutate(Stage = case_when(
    Stage %in% c("Stage 1", "Stage 2", "Stage 3", "Stage 4") ~ Stage,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Stage))

ClinicalFeatures <- ClinicalFeatures %>%
  filter(race != "not reported")

merged_data <- ClinicalFeatures[,Prognostic_gene]%>% t() %>%as.data.frame()

head(merged_data) 
annotation <- ClinicalFeatures %>%
  select( RiskGroup, gender, T.stage,N.stage, Stage, race, Age)
my_colors <- colorRampPalette(c("blue", "white", "red"))(50)

pdf('./result/figure/07_Fig7D_Express_heatmap.pdf',width = 9,height = 12,onefile = F)
pheatmap(merged_data, 
         annotation_col = annotation, 
         color = my_colors,
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         scale = "row")
dev.off()

ClinicalFeatures$OS <-ifelse(ClinicalFeatures$OS == "Alive",0,1)

saveRDS(ClinicalFeatures,file = './data/07_ClinicalFeatures.rds')


ClinicalFeatures <-readRDS('./data/07_ClinicalFeatures.rds')


head(ClinicalFeatures)


############################################################
#
#-------------------生存图。——————————
#
############################################################

ClinicalFeatures$AgeGroup <- ifelse(ClinicalFeatures$A1_Age > 65, ">65", "≤65")
ClinicalFeatures$TumorGrade <- ifelse(ClinicalFeatures$ajcc_pathologic_stage %in% c("Stage I", 'Stage IA','Stage IB',"Stage IIA","Stage IIB"), "I-II", "III-IV")
ClinicalFeatures$NStage <- ifelse(ClinicalFeatures$N.stage == "N0", "N0", "N1-3")
ClinicalFeatures$TStage <- ifelse(ClinicalFeatures$T.stage %in% c("T1", "T2"), "T1-2", "T3-4")
ClinicalFeatures$PathologicStage <- ifelse(ClinicalFeatures$Stage %in% c("Stage 1", "Stage 2"), "I-II", "III-IV")




# 确保时间变量是数值型
ClinicalFeatures$OS.time <- as.numeric(ClinicalFeatures$OS.time)

# 创建生存对象
surv_object <- Surv(time = ClinicalFeatures$OS.time, event = ClinicalFeatures$OS)

# 查看生存对象
print(surv_object)


# 根据不同亚组进行生存分析
fit_age <- survfit(surv_object ~ AgeGroup, data = ClinicalFeatures)
fit_gender <- survfit(surv_object ~ gender, data = ClinicalFeatures)
fit_tumor_grade <- survfit(surv_object ~ TumorGrade, data = ClinicalFeatures)
fit_n_stage <- survfit(surv_object ~ NStage, data = ClinicalFeatures)
fit_t_stage <- survfit(surv_object ~ TStage, data = ClinicalFeatures)
fit_pathologic_stage <- survfit(surv_object ~ PathologicStage, data = ClinicalFeatures)

# 创建多个生存曲线图，并设置颜色
plot_age <- ggsurvplot(fit_age, data = ClinicalFeatures, pval = TRUE, title = "Survival by Age Group",
                       palette = c("#EDAA25", "#C43302"))
plot_gender <- ggsurvplot(fit_gender, data = ClinicalFeatures, pval = TRUE, title = "Survival by Gender",
                          palette = c("#EDAA25", "#C43302"))
plot_tumor_grade <- ggsurvplot(fit_tumor_grade, data = ClinicalFeatures, pval = TRUE, title = "Survival by Tumor Grade",
                               palette = c("#EDAA25", "#C43302"))
plot_n_stage <- ggsurvplot(fit_n_stage, data = ClinicalFeatures, pval = TRUE, title = "Survival by N Stage",
                           palette = c("#EDAA25", "#C43302"))
plot_t_stage <- ggsurvplot(fit_t_stage, data = ClinicalFeatures, pval = TRUE, title = "Survival by T Stage",
                           palette = c("#EDAA25", "#C43302"))
plot_pathologic_stage <- ggsurvplot(fit_pathologic_stage, data = ClinicalFeatures, pval = TRUE, title = "Survival by Pathologic Stage",
                                    palette = c("#EDAA25", "#C43302"))

# 提取ggplot对象
p1 <- plot_age$plot
p2 <- plot_gender$plot
p3 <- plot_tumor_grade$plot
p4 <- plot_n_stage$plot
p5 <- plot_t_stage$plot
p6 <- plot_pathologic_stage$plot

# 使用cowplot组合图形
combined_plot <- plot_grid(p1, p2, p3, p4, p5, p6, labels = "AUTO", ncol = 2)

combined_plot
# 保存组合图形
ggsave("./result/figure/07_Fig7A_combined_survival_plots.pdf", combined_plot, width = 12, height = 18)















 

