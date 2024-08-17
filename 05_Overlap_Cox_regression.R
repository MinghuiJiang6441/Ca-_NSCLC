setwd("/data/home/jiangminghui/NSCLC")
library(glmnet)
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
library(VennDiagram)
#install.packages("venn")
library(ggvenn)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(devtools)
#install_github('ievaKer/aPEAR')
library(DOSE)
library(forestplot)
library(forestploter)
library(gridExtra)
#library(aPEAR)
library(survminer)
library(ggrisk)
library(survival)
############################################################
#
#--------------------单因素cox————————————————————————————
#
############################################################

#table(Overlap_genes %in%rownames(expr_combat))

expr <- readRDS('./data/02_TCGA_expr.rds')

Overlap_genes <-readRDS('./result/data/Overlap_genes.rds')

TCGA_NSCLC_clin <- readRDS('./data/02_TCGA_NSCLC_clin.rds')

cox_exp <- expr[Overlap_genes,]

cox_exp <- log2(cox_exp + 1)

range(cox_exp)

cox_data <-cbind(TCGA_NSCLC_clin,t(cox_exp))

# 将"Alive"和"Dead"转换为0和1
cox_data$A7_OS.event <- ifelse(cox_data$A7_OS.event == "Alive", 0, 1)

cox_data <- cox_data[-which(is.na(cox_data$A7_OS.time )),]
cox_data$A7_OS.time <- as.numeric(cox_data$A7_OS.time)

cox_data <- cox_data [-which(cox_data$A7_OS.time ==0) ,]
# 存储结果的列表
results <- list()

cox_data <-readRDS('./result/data/05_cox_data.rds')
# 单因素Cox回归分析
for (gene in Overlap_genes) {
  formula <- as.formula(paste("Surv(A7_OS.time, A7_OS.event) ~", gene))
  cox_model <- coxph(formula, data = cox_data)
  summary_cox <- summary(cox_model)
  p_value <- summary_cox$coefficients[5]
  
  if (p_value < 0.01) {
    results[[gene]] <- summary_cox
  }
}
forest_data <- do.call(rbind, lapply(results, function(x) {
  coef <- x$coefficients[1]
  se <- x$coefficients[3]
  lower <- coef - 1.96 * se
  upper <- coef + 1.96 * se
  c(coef, lower, upper)
}))


dim(forest_data)

forest_data <- as.data.frame(forest_data)

head(forest_data)
rownames(forest_data) <- names(results)

forest_data$Gene <- rownames(forest_data)

rownames(forest_data) <- names(results)

head(forest_data)

forestplot(labeltext = forest_data$Gene,
           mean = forest_data$HR,
           lower = forest_data$Lower,
           upper = forest_data$Upper,
           xlab = "Hazard Ratio",
           title = "Forest Plot of Hazard Ratios and p-values",
           new_page = TRUE,
           is.summary = FALSE,
           boxsize = 0.2,
           col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"))







head(cox_data)
univ_formulas <- sapply(Overlap_genes,
                        function(x) as.formula(paste('Surv(A7_OS.time, A7_OS.event)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data =  cox_data )})      

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR,signif(x$coef[2], digits=2),HR.confint.lower,HR.confint.upper)
                         names(res)<-c("p.value","HR (95% CI for HR)",'HR','lower','upper')
                         return(res)
                       }) 

univ_results <- t(as.data.frame(univ_results, check.names = FALSE))
univ_results <-as.data.frame(univ_results,stringsAsFactors=F) 
univ_results$HR <- as.numeric(univ_results$HR)
univ_results$lower <- as.numeric(univ_results$lower)
univ_results$upper <- as.numeric(univ_results$upper)
tabletext<-cbind(c(NA,"module",rownames(univ_results)),
                 c(NA,"Hazard Ratio(95% CI)",univ_results$`HR (95% CI for HR)`),
                 c(NA,"pValue",ifelse(univ_results$p.value<0.001,"P < 0.001",univ_results$p.value))
)

cochrane_from_rmeta <- 
  structure(list(
    mean  = univ_results$HR, 
    lower = univ_results$lower,
    upper = univ_results$upper),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -10L), 
    class = "data.frame")

colnames(tabletext) <-tabletext[2,]
tabletext <-tabletext[-c(1:2) ,]

save(tabletext,tabletext,cochrane_from_rmeta,file = './result/data/05_cox_forestplot.RData')

head(tabletext)

pdf('./result/figure/05_Fig5A_Cox_forestPlot.pdf',onefile = F)
forestplot(tabletext, 
              cochrane_from_rmeta,new_page = TRUE,
              #is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
             # hrzl_lines = list("3" = gpar(lty=2), 
              #                  "11" = gpar(lwd=1, columns=1:4, col = "#000044")),
              clip=c(0.1,2.5),
              title=" Univariate OS Cox Forest plot",
              
              xlog=TRUE,
              col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
dev.off()





genes <- c('TMEM132A','SLC16A1','ARNTL2','CTSV')
# 循环处理每个基因

  # 拟合Cox回归模型
  fit.cox <- coxph(as.formula(paste('Surv(A7_OS.time, A7_OS.event) ~', 'CTSV')), data = cox_data)
  
  # 提取模型结果
  tidy_fit <- broom::tidy(fit.cox, exponentiate = TRUE, conf.int = TRUE)
  
  # 进行比例风险假设检验
  ftest <- cox.zph(fit.cox)
  
  # 生成图形
  plot <- ggcoxzph(ftest)
  
  # 将图形添加到列表中
  pdf(paste0('./result/figure/05_Fig5B_','CTSV','.pdf'),onefile = FALSE,width = 6,height = 8)
  
 plot
 dev.off()
 
 saveRDS(cox_data,file = './result/data/05_cox_data.rds')
 
 #Prognostic_gene <- c('TMEM132A','SLC16A1','ARNTL2','CTSV')



