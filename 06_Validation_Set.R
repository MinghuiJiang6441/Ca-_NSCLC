############################################################
#
#--------------------GSE19188———————————————————————————
#
############################################################

bothsample <- GSE19188_clin[which(GSE19188_clin$tissue.type =='tumor'),] %>%rownames()
head(bothsample)

GSE19188_exp <-GSE19188_exp[,bothsample]
GSE19188_clin <- GSE19188_clin[bothsample,]


x <- GSE19188_exp[Prognostic_gene,] %>%t()

GSE19188_Lasso_data <- cbind(GSE19188_clin[,c(11,12)],x)

GSE19188_Lasso_data <- subset(GSE19188_Lasso_data,GSE19188_Lasso_data$overall.survival !='Not available')

GSE19188_Lasso_data$overall.survival <- as.double(GSE19188_Lasso_data$overall.survival )
GSE19188_Lasso_data$status <- ifelse(GSE19188_Lasso_data$status == "alive", 0, 1)


GSE19188_Lasso_data$status <- as.double(GSE19188_Lasso_data$status )


coef_selectGene

head(GSE19188_Lasso_data)

coef_selectGene <- as.vector(coef_selectGene)

GSE19188_Lasso_data$RiskScore <- rowSums(GSE19188_Lasso_data[, -c(1,2)] * coef_selectGene)

threshold <- median(GSE19188_Lasso_data$RiskScore)

GSE19188_Lasso_data$RiskGroup <- ifelse(GSE19188_Lasso_data$RiskScore > threshold, "High_Risk", "Low_Risk")

pca_prcp <-  GSE19188_Lasso_data[,3:13]  %>% prcomp 


pca_prcp$x %>% .[,1:2]  %>% head
pca_prcp_contrib <- pca_prcp$sdev %>% .^2 %>% {./sum(.) * 100} %>% .[1:2] %>% signif(digits = 4)
pca_prcp_contrib

# 确保pca_prcp$x存在并且包含前两个主成分
pca_data <- pca_prcp$x[, 1:2] %>% as.data.frame()

# 确保risk_data$RiskGroup存在并且长度匹配
pca_data <- cbind(pca_data, Groups = GSE19188_Lasso_data$RiskGroup)

pdf('./result/figure/06_fig6E_Risk_PCA_GSE19188.pdf',width = 5,height = 5)

ggplot(pca_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Groups)) +
  labs(title = "PCA of NSCLC Samples",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_color_manual(values = c("Low_Risk" = "#0455BF", "High_Risk" = "#8C3414")) +
  
  theme_minimal()
dev.off()

GSE19188_Lasso_data$overall.survival <- GSE19188_Lasso_data$overall.survival*12

cutoff<-surv_cutpoint(GSE19188_Lasso_data, #数据集
                      time="overall.survival",#“ ”里写数据集时间变量
                      event="status",##“ ”里数据集结局变量名称
                      variables='RiskScore'
);

summary(cutoff) #输出结果
groups<-surv_categorize(cutoff)

table(groups$RiskScore)



fit <- survfit(Surv(overall.survival , status) ~ RiskGroup, data = GSE19188_Lasso_data)

pdf('./result/figure/06_Fig6C_KM_curve.pdf',width = 5,height = 7,onefile = FALSE)
ggsurvplot(
  fit, 
  data = GSE19188_Lasso_data, 
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


with(GSE19188_Lasso_data,
     ROC_riskscore <<- timeROC(T = overall.survival,
                               delta = status,
                               marker = RiskScore,
                               cause = 1,
                               weighting = "marginal",
                               times = c(365,1080,1800),
                               ROC = TRUE,
                               iid = TRUE)
)

plot(ROC_riskscore, time = 365, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 1080, col = "blue", add = T)
plot(ROC_riskscore, time = 1800, col = "purple", add = T)
legend("bottomright",c("1-Year","3-Year","5-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
text(0.5,0.15,paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
text(0.5,0.1,paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3)))


############################################################
#
#--------------------GSE19188———————————————————————————
#
############################################################

dim(GSE37745_exp)
head(GSE37745_clin)

a <- GSE37745_exp[Prognostic_gene,] %>% t()

colnames(GSE37745_clin)

GSE37745_data <-cbind(GSE37745_clin[,c(17,18)],a)

names(GSE37745_data)[c(1,2)] <- c('OS.time','OS')

head(GSE37745_data,2)

GSE37745_data$OS.time <- as.double(GSE37745_data$OS.time )

GSE37745_data$OS <- ifelse(GSE37745_data$OS == "no", 0, 1)

GSE37745_data$OS <- as.double(GSE37745_data$OS )


GSE37745_data$RiskScore <- rowSums(GSE37745_data[, -c(1,2)] * coef_selectGene)

threshold <- median(GSE19188_Lasso_data$RiskScore)

GSE37745_data$RiskGroup <- ifelse(GSE37745_data$RiskScore > threshold, "High_Risk", "Low_Risk")

pca_prcp <-  GSE37745_data[,3:13] %>% scale() %>% prcomp 


pca_prcp$x %>% .[,1:2]  %>% head
pca_prcp_contrib <- pca_prcp$sdev %>% .^2 %>% {./sum(.) * 100} %>% .[1:2] %>% signif(digits = 4)
pca_prcp_contrib

# 确保pca_prcp$x存在并且包含前两个主成分
pca_data <- pca_prcp$x[, 1:2] %>% as.data.frame()

# 确保risk_data$RiskGroup存在并且长度匹配
pca_data <- cbind(pca_data, Groups = GSE37745_data$RiskGroup)

pdf('./result/figure/06_Fig6E_PCA_GSE37745.pdf',width = 5,height = 5,onefile = FALSE)
ggplot(pca_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Groups)) +
  labs(title = "PCA of NSCLC Samples",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_color_manual(values = c("Low_Risk" = "#0455BF", "High_Risk" = "#8C3414")) +
  
  theme_minimal()
dev.off()

fit <- survfit(Surv(OS.time , OS) ~ RiskGroup, data = GSE37745_data)

pdf('./result/figure/06_Fig6G_KM_curve_GSE37745.pdf',width = 5,height = 7,onefile = FALSE)

ggsurvplot(
  fit, 
  data = GSE37745_data, 
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

formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(Prognostic_gene, collapse = "+")))

# 拟合Cox回归模型
res.cox <- coxph(formula, data = GSE37745_data)

res.cox

pdf('./result/figure/06_Fig6F_Risk_score_GSE37745.pdf',width = 8,height = 15,onefile = FALSE)

ggrisk(res.cox,
       cutoff.value='median',
       cutoff.x = 145,
       cutoff.y = -0.8,
       code.0 = 'Still Alive',
       code.1 = 'Already Dead',
       code.highrisk = 'High Risk',
       code.lowrisk = 'Low Risk',
       title.A.ylab='Risk Score',
       title.B.ylab='Survival Time(year)',
       title.A.legend='Risk Group',
       title.B.legend='Status',
       title.C.legend='Expression',
       color.A=c(low ='#0455BF',high='#8C3414'),#A图中点的颜色
       color.B=c(code.0='#049DD9',code.1='#7D73BF'), #B图中点的颜色
       #  color.C=c(low='#BEB8C6',median='#735F95',high='#27204C'), #C图中热图颜色,
       heatmap.genes=Prognostic_gene
)

dev.off()

with(GSE37745_data,
     ROC_riskscore <<- timeROC(T = OS.time,
                               delta = OS,
                               marker = RiskScore,
                               cause = 1,
                               weighting = "marginal",
                               times = c(365,1080,1800),
                               ROC = TRUE,
                               iid = TRUE)
)
pdf('./result/figure/06_Fig6H_ROc_curve_GSE37745.pdf',width = 6,height = 6,onefile = FALSE)

plot(ROC_riskscore, time = 365, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 1080, col = "blue", add = T)
plot(ROC_riskscore, time = 1800, col = "purple", add = T)
legend("bottomright",c("1-Year","3-Year","5-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
text(0.5,0.15,paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
text(0.5,0.1,paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3)))

dev.off()
