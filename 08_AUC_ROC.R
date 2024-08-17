library(MASS)
library(survivalROC)
library(compareC)
library("doParallel")
remotes::install_github("nanxstats/hdnom")
library("hdnom")
library(rms)
library(gridExtra)
library(pROC)
library(forestplot)
library(survminer)


ClinicalFeatures$OS.time <- as.double(ClinicalFeatures$OS.time )
ClinicalFeatures$OS <- as.double(ClinicalFeatures$OS )

colnames(ClinicalFeatures)

head(ClinicalFeatures,2)
############################################################
#
#-------------------单因素——————————————————————————
#
############################################################
CoxData  <- ClinicalFeatures %>%
  mutate_at(vars(A1_Age, RiskScore, race, gender, Stage, N.stage, T.stage), ~ as.numeric(as.factor(.)))
#假设我们要对如下5个特征做单因素cox回归分析
covariates <- c('A1_Age','RiskScore',"race", "gender",  "Stage", "N.stage", "T.stage")
#分别对每一个变量，构建生存分析的公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time, OS)~', x)))

#循环对每一个特征做cox回归分析
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = CoxData)})

#res.cox <- coxph(Surv(OS.time, OS ) ~ RiskScore, data = ClinicalFeatures)

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

#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

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

pdf('./result/figure/08_Fig8A_singlevariate_cox.pdf',onefile = F)

forestplot(tabletext,cochrane_from_rmeta[1:7,],xlog=TRUE)
dev.off()

############################################################
#
#-------------------多因素———————————————————————————
#
############################################################
bigmodel <-
  coxph(Surv(OS.time, OS) ~ A1_Age + RiskScore+race+gender+Stage+ N.stage+T.stage,
        data = CoxData )

pdf('./result/figure/08_Fig8A_Multivariate_cox.pdf',onefile = F)
ggforest(bigmodel)
dev.off()



stepwise_model <- stepAIC(bigmodel, direction = "both")

summary(stepwise_model) #独立预后因子

str(stepwise_model)

############################################################
#
#-------------------计算预后模型的AUC值——————————————————————
#
############################################################
# 预测风险评分
var <- c('A1_Age','RiskScore','gender','Stage','T.stage')
x <- as.matrix(CoxData[,var])
time <- CoxData$OS.time
event <- CoxData$OS
y <- survival::Surv(time, event)
registerDoParallel(detectCores())

fit <- fit_aenet(x, y, nfolds = 10, rule = "lambda.1se", seed = c(5, 7), parallel = TRUE)


fit <- fit_aenet(x, y, nfolds = 10, rule = "lambda.1se", seed = c(5, 7), parallel = TRUE)










############################################################
#
#------------------柱状图——————————————————————
#
############################################################
# 预测风险评分
multiv_cox <- coxph(Surv(OS.time, OS) ~ A1_Age + RiskScore+gender+Stage+T.stage, data = CoxData)

risk_score <- predict(multiv_cox, type = "risk")

# 计算模型的一致性指数（C-index）
c_index_model <- concordance.index(risk_score, data$time, data$status)$c.index



dd <- datadist(CoxData)
options(datadist = "dd")

mul_cox_2 <- cph(Surv(OS.time, OS) ~ A1_Age + RiskScore+gender+Stage+T.stage, data = CoxData, x = TRUE, y = TRUE, surv = TRUE)

sur <- Survival(mul_cox_2)

sur1 <- function(x) sur(365, x)   # 1年生存
sur2 <- function(x) sur(730, x)   # 2年生存
sur3 <- function(x) sur(1095, x)  # 3年生存

nom <- nomogram(mul_cox_2, 
                fun = list(sur1, sur2, sur3),  # 为每个时间点指定生存函数
                fun.at = c(0.05, seq(0.1, 0.9, by = 0.05), 0.95),  # 设置诺莫图上的生存概率点
                funlabel = c("1 year survival", "2 years survival", "3 years survival"))  # 设置生存时间标签

pdf('./result/figure/08_Fig8D_nomograms_Linechart.pdf',width = 12,height = 6)
plot(nom)
dev.off()

cal1 <- calibrate(mul_cox_2, 
                  cmethod = 'KM',   # 表示使用Kaplan-Meier（KM）估计方法进行校准
                  method = "boot",  # 表示使用自助法（Bootstrap）进行校准，Bootstrap 是一种统计方法，它通过从原始数据中有放回地进行重采样来估计参数的不确定性和分布。在这里，Bootstrap 用于生成多个随机样本来估计校准曲线的分布，以便获得更可靠的校准结果。
                  u = 365,          # 设置时间间隔，需要与之前模型中定义的time.inc一致
                  m = 66,           # 每次抽样的样本量，根据样本量来确定，标准曲线一般将所有样本分为3组（在图中显示3个点）
                  B = 1000)         # 抽样次数

par(mar = c(6, 6, 3, 3))
 plot(cal1,                          # 绘制校准曲线的数据
     lwd=1,                         # 线条宽度为1
     lty=1,                         # 线条类型为1（实线）
     conf.int=T,                    # 是否显示置信区间
     errbar.col="blue3",            # 直线和曲线的误差线颜色设置为蓝色
     col="red3",                    # 校准曲线的颜色设置为红色
     xlim=c(0,1),                   # x轴的限制范围，从0到1
     ylim=c(0,1),                   # y轴的限制范围，从0到1
     xlab="Nomogram-Predicted Probability of 1-Year OS",  # x轴标签
     ylab="Actual 1-Year DFS (proportion)",                # y轴标签
     subtitles = F)                 # 不显示副标题



 # Function to calibrate and plot
 calibrate_and_plot <- function(model, time, xlab, ylab) {
   cal <- calibrate(model, cmethod = 'KM', method = "boot", u = time, m = 66, B = 1000)
   plot(cal, lwd=1, lty=1, conf.int=T, errbar.col="blue3", col="red3", xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, subtitles = F)
 }
 
 # 生成图形对象
 p1 <- recordPlot(calibrate_and_plot(mul_cox_2, 365, "Nomogram-Predicted Probability of 1-Year OS", "Actual 1-Year DFS (proportion)"))
 p2 <- recordPlot(calibrate_and_plot(mul_cox_2, 365*3, "Nomogram-Predicted Probability of 3-Year OS", "Actual 3-Year DFS (proportion)"))
 p3 <- recordPlot(calibrate_and_plot(mul_cox_2, 365*5, "Nomogram-Predicted Probability of 5-Year OS", "Actual 5-Year DFS (proportion)"))

 grid.arrange(p1, p2, p3, ncol=3)
 

 # 创建一个1行3列的布局
 pdf('./result/figure/08_Fig8E_nCalibration_curve.pdf',width = 9,height = 3)
 layout(matrix(c(1, 2, 3), 1, 3, byrow = TRUE))
 
 p1
p2
p3
dev.off()

with(CoxData,
     ROC_riskscore <<- timeROC(T = OS.time,
                               delta = OS,
                               marker = RiskScore,
                               cause = 1,
                               weighting = "marginal",
                               times = c(365,1080,1800),
                               ROC = TRUE,
                               iid = TRUE)
)


head(ROC_riskscore)

op <- par(no.readonly = TRUE)  # 保存当前参数
# 绘图代码
par(op)  # 恢复原始参数

graphics.off()

pdf('./result/figure/08_Fig8F.pdf',width = 5,height = 5,onefile = F)

plot.roc(CoxData$OS, CoxData$RiskScore, # data
         
         percent=TRUE, # show all values in percent
         
         partial.auc=c(100, 90), partial.auc.correct=TRUE, # define a partial AUC (pAUC)
         
         print.auc=TRUE, #display pAUC value on the plot with following options:
         
         print.auc.pattern="Corrected pAUC (100-90%% SP):\n%.1f%%", print.auc.col="#1c61b6",
         
         auc.polygon=TRUE, auc.polygon.col="#1c61b6", # show pAUC as a polygon
         
         max.auc.polygon=TRUE, max.auc.polygon.col="#1c61b622", # also show the 100% polygon
         
         main="Partial AUC (pAUC)")

plot.roc(CoxData$OS, CoxData$RiskScore,
         
         percent=TRUE, add=TRUE, type="n", # add to plot, but don't re-add the ROC itself (useless)
         
         partial.auc=c(100, 90), partial.auc.correct=TRUE,
         
         partial.auc.focus="se", # focus pAUC on the sensitivity
         
         print.auc=TRUE, print.auc.pattern="Corrected pAUC (100-90%% SE):\n%.1f%%", print.auc.col="#008600",
         
         print.auc.y=40, # do not print auc over the previous one
         
         auc.polygon=TRUE, auc.polygon.col="#008600",
         
         max.auc.polygon=TRUE, max.auc.polygon.col="#00860022")

dev.off()








