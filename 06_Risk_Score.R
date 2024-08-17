library(magrittr)
library(ggplot2)
library(ggrisk)
setwd("/data/home/jiangminghui/NSCLC")
library(timeROC)
library(survival)


TCGA_NSCLC_clin <- readRDS('./data/02_TCGA_NSCLC_clin.rds')

dim(TCGA_NSCLC_clin)

coef_selectGene <- coef(cv.fit, s="lambda.min")[select.varialbes,]

exp_subset <- expr[Prognostic_gene,] %>% t()

risk_data <- cbind(TCGA_NSCLC_clin[,c(4,5)],exp_subset)

#Risk_score <- apply(exp_subset,2,function(x){
#  sum(x*coef_subset)
#})

#TCGA_NSCLC_clin$Risk_score <- Risk_score



#saveRDS(TCGA_NSCLC_clin,'./data/02_TCGA_NSCLC_clin.rds')

head(risk_data)

risk_data$A7_OS.event <- ifelse(risk_data$A7_OS.event == "Alive", 0, 1)

risk_data$A7_OS.time <- as.numeric(risk_data$A7_OS.time)

risk_data <- risk_data[-which(is.na(risk_data$A7_OS.time )),]

risk_data <- risk_data [-which(risk_data$A7_OS.time ==0) ,]

formula <- as.formula(paste("Surv(A7_OS.time, A7_OS.event) ~", paste(Prognostic_gene, collapse = "+")))

# 拟合Cox回归模型
res.cox <- coxph(formula, data = risk_data)

res.cox

saveRDS(res.cox,file = './result/data/06_Fig6B_res.cox.rds')
#install.packages('ggrisk')
pdf('./result/figure/06_Fig6B_Risk_score.pdf',width = 8,height = 15,onefile = FALSE)
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

str(res.cox)


a <- ggrisk(res.cox,
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
            heatmap.genes=Prognostic_gene)
str(res.cox)


ggrisk

colnames(risk_data)


risk_data$RiskScore <- predict(res.cox, type = "lp")
risk_data$RiskScore <- rowSums(risk_data[, -c(1,2,14,15)] * coef_selectGene)


threshold <- median(risk_data$RiskScore)

risk_data$RiskGroup <- ifelse(risk_data$RiskScore > threshold, "High_Risk", "Low_Risk")

head(risk_data)

#features <- risk_data[, -c(1, which(names(risk_data) %in% c("risk_score", "risk_group")))]

saveRDS(risk_data,'./result/data/risk_data.rds')
pca_prcp <-  risk_data[,3:13] %>% scale() %>% prcomp 


pca_prcp$x %>% .[,1:2]  %>% head
pca_prcp_contrib <- pca_prcp$sdev %>% .^2 %>% {./sum(.) * 100} %>% .[1:2] %>% signif(digits = 4)
pca_prcp_contrib

# 确保pca_prcp$x存在并且包含前两个主成分
pca_data <- pca_prcp$x[, 1:2] %>% as.data.frame()

# 确保risk_data$RiskGroup存在并且长度匹配
pca_data <- cbind(pca_data, Groups = risk_data$RiskGroup)



pdf('./result/figure/06_fig6A_Risk_PCA.pdf',width = 5,height = 5)
ggplot(pca_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Groups)) +
  labs(title = "PCA of NSCLC Samples",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_color_manual(values = c("Low_Risk" = "#0455BF", "High_Risk" = "#8C3414")) +
  
  theme_minimal()
dev.off()


############################################################
#
#-------------------km曲线———————————————————————————
#
############################################################

head(risk_data,2)

fit <- survfit(Surv(A7_OS.time , A7_OS.event) ~ RiskGroup, data = risk_data)

pdf('./result/figure/06_Fig6C_KM_curve.pdf',width = 5,height = 7,onefile = FALSE)
ggsurvplot(
  fit, 
  data = risk_data, 
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

############################################################
#
#-------------------ROC曲线———————————————————————————
#
############################################################

with(risk_data,
     ROC_riskscore <<- timeROC(T = A7_OS.time,
                               delta = A7_OS.event,
                               marker = RiskScore,
                               cause = 1,
                               weighting = "marginal",
                               times = c(365,1080,1800),
                               ROC = TRUE,
                               iid = TRUE)
)

pdf('./result/figure/06_Fig6D_ROc_curve.pdf',width = 6,height = 6,onefile = FALSE)
plot(ROC_riskscore, time = 365, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 1080, col = "blue", add = T)
plot(ROC_riskscore, time = 1800, col = "purple", add = T)
legend("bottomright",c("1-Year","3-Year","5-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
text(0.5,0.15,paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
text(0.5,0.1,paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3)))

dev.off()




















