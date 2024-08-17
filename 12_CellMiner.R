library(pRRophetic)
library(ggplot2)
library(ggpubr)
library(parallel)

data(bortezomibData)
dim(exprDataBortezomib)
table(studyResponse)
predictedPtype <- pRRopheticPredict(testMatrix=as.matrix(expr[,CoxData$Sample]), 
                                    drug="GDC.0449",
                                    tissueType = "all", 
                                    batchCorrect = "eb",
                                    selection=1,
                                    dataset = "cgp2014")

#将数据类型转化为数据框
v<- as.data.frame(predictedPtype)
#将行名转化为列名
predictedPtype <- rownames_to_column(predictedPtype, var = "Sample")

predictedPtype$group <-CoxData$RiskGroup

ggplot(data=predictedPtype,aes(x=group,y= predictedPtype))+
  geom_violin(aes(fill=group))+
  geom_boxplot(width=0.2)+
  # 显著差异分析
  geom_signif(comparisons=list(c("High_Risk", "Low_Risk")),
              test=wilcox.test,
              step_increase=0.2,
              map_signif_level=TRUE) +
  scale_fill_npg()+
  labs(x="",y="Bortezomib IC50")+
  theme_classic() 

data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016) # drugData2016
data(cgp2016ExprRma) # cgp2016ExprRma
data(bortezomibData)
#提取cgp2016的所有药物名字

possibleDrugs2016 <- unique( drugData2016$Drug.name)

plot_drug_response <- function(expr, CoxData, drug_name) {
  predictedPtype <- pRRopheticPredict(testMatrix=as.matrix(expr[,CoxData$Sample]), 
                                      drug=drug_name,
                                      tissueType = "all", 
                                      batchCorrect = "eb",
                                      selection=1,
                                      dataset = "cgp2014")
  
  # 将数据类型转化为数据框
  predictedPtype <- as.data.frame(predictedPtype)
  # 将行名转化为列名
  predictedPtype <- rownames_to_column(predictedPtype, var = "Sample")
  
  predictedPtype$group <- CoxData$RiskGroup
  
  # 绘制图像
  p <- ggplot(data=predictedPtype, aes(x=group, y=predictedPtype)) +
    geom_violin(aes(fill=group)) +
    geom_boxplot(width=0.2) +
    # 显著差异分析
    geom_signif(comparisons=list(c("High_Risk", "Low_Risk")),
                test=wilcox.test,
                step_increase=0.2,
                map_signif_level=TRUE) +
    scale_fill_npg() +
    labs(x="", y=paste(drug_name, "IC50")) +
    theme_classic()
  

  pdf(paste0('./result/figure/12_Fig12A_',drug_name,'_boxplot.pdf'), onefile = FALSE)
  print(p)
  dev.off()
}

# 使用函数
plot_drug_response(expr, CoxData,  "PLX4720"  )
                           
drug_list <- c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", 
               "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", 
               "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", 
               "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", 
               "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", 
               "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", 
               "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", 
               "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", 
               "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", 
               "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", 
               "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", 
               "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", 
               "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", 
               "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", 
               "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", 
               "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", 
               "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", 
               "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", 
               "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", 
               "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", 
               "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")

  length(drug_list)                        
  # 随机选择10种药物
  set.seed(123) # 设置随机种子以确保结果可重复
  selected_drugs <- sample(drug_list, 10)
  
  # 初始化空数据框
  all_predictions <- data.frame(Sample = CoxData$Sample)
  
  # 循环预测药物敏感性并合并数据框
  for (drug in selected_drugs) {
    predictedPtype <- pRRopheticPredict(testMatrix=as.matrix(expr[,CoxData$Sample]), 
                                        drug=drug,
                                        tissueType = "all", 
                                        batchCorrect = "eb",
                                        selection=1,
                                        dataset = "cgp2014")
    predictedPtype <- as.data.frame(predictedPtype)
    predictedPtype <- rownames_to_column(predictedPtype, var = "Sample")
    
    all_predictions <- merge(all_predictions, predictedPtype, by = "Sample", all = TRUE)
  }
  
  selected_drugs                     
colnames(all_predictions) <-c('Sample',selected_drugs)
  head(all_predictions)
  
  head(CoxData)
  
  # 计算Spearman相关系数
  correlation_matrix <- corr.test(all_predictions[,-1], CoxData[,c(Prognostic_gene,'RiskScore')], method = "spearman")

  my_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
  
  # 绘制热图并改变颜色
  pdf('./result/figure/12_Fig12B_drug_sensitivity_correlation_heatmap.pdf', onefile = FALSE)
  pheatmap(correlation_matrix$r, 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           display_numbers = matrix(ifelse(correlation_matrix$p < 0.05, "*", ""), nrow(correlation_matrix$p)), 
           main = "Spearman Correlation between Drug Sensitivity and Gene Expression",
           color = my_colors)
  dev.off()