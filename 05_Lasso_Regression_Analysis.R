library(glmnet)
library(ggradar)
library(scales)
library(tibble)
library(fmsb)


############################################################
#
#--------------------Lasso———————————————————————————
#
############################################################
head(Lasso_data )

#Lasso_data <- cox_data[,c('A7_OS.time','A7_OS.event',Prognostic_gene)]
which(colnames(cox_data) =='NAPSA')

dim(cox_data)

Lasso_data <-cox_data[,c(4,5,c(11:56))]

Lasso_data$A7_OS.time <- as.double(Lasso_data$A7_OS.time )
Lasso_data$A7_OS.event <- as.double(Lasso_data$A7_OS.event )
#Lasso_data <- Lasso_data[-which(Lasso_data$A7_OS.time==0),]

x <- Lasso_data[,-c(1,2)] %>%scale() %>%as.data.frame()


y <- data.matrix(Surv(Lasso_data$A7_OS.time,Lasso_data$A7_OS.event))


fit <-glmnet(x,y,family = "cox",alpha = 1)

pdf('./result/figure/05_Fig5C_Lasso_path_diagram.pdf',width = 5,height = 5,onefile = FALSE)
plot(fit,xvar="lambda",label=T)
dev.off()



set.seed(4)#不要改！！！！！！！
# 使用cv.glmnet函数
cv.fit <- cv.glmnet(as.matrix(x), y, family = 'cox', alpha = 1, nfolds = 10)
best_lambda <- cv.fit$lambda.min
pdf('./result/figure/05_Fig5C_Lasso_path_diagram2.pdf',onefile = FALSE)

plot(cv.fit,xvar='lambda',label=TRUE)
abline(v=log(best_lambda),col='red',lty=2)
dev.off()

select.varialbes = rownames(as.data.frame(which(coef(cv.fit, s = "lambda.min")[,1]!=0)))

Prognostic_gene <- select.varialbes

length(Prognostic_gene)

saveRDS(cv.fit,file = './result/data/05_cv.fit.rds')

saveRDS(Prognostic_gene,file = './result/data/05_Prognostic_gene.rds') 

############################################################
#
#--------------------雷达图————————————————————————————
#
############################################################
coef_selectGene <- coef(cv.fit, s="lambda.min")[Prognostic_gene,] 

#%>%t()%>%as.data.frame()


ggradar(coef_selectGene) 

head(mtcars)

coef_selectGene <- data.frame(
  TMPRSS2 = -0.010803391,
  GGTLC1 = -0.009142888,
  AGR3 = 0.008908738,
  ALDH3B1 = 0.152304168,
  SGMS2 = 0.107053401,
  EPB41L5 = -0.024429694,
  UNC13B = -0.009730462,
  TMEM163 = -0.088898294,
  GPD1L = -0.170060373,
  SLC16A1 = 0.040558301,
  ARNTL2 = 0.045297357
)

coef_selectGene <- rbind(rep(max(coef_selectGene), 11), rep(min(coef_selectGene), 11), coef_selectGene)

pdf('./result/figure/05_Fig5D_Coefficients_redar_plot.pdf',width = 6,height = 6,onefile = FALSE)
# 绘制雷达图
radarchart(coef_selectGene, axistype = 1,
           pcol = rgb(0.2, 0.5, 0.5, 0.9), pfcol = rgb(0.2, 0.5, 0.5, 0.5), plwd = 2,
           cglcol = "grey", cglty = 1, axislabcol = "grey", caxislabels = seq(-0.2, 0.2, 0.1), cglwd = 0.8,
           vlcex = 0.8)

dev.off()

#验证集


save.image('./0809.RData')
q()










