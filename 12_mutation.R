library(easyTCGA)
library(maftools)
library(dplyr)

CoxData <- readRDS(file = './result/data/CoxData.rds')
# 直接读取
getsnvmaf("TCGA-LUAD")
getsnvmaf("TCGA-LUSC")

load('./output_snv/TCGA-LUAD_matched_maf_and_clin.rdata')
load('./output_snv/TCGA-LUSC_matched_maf_and_clin.rdata')
#LUAD_maf <- read.maf(data,clin_snv,isTCGA = T)

LUAD_data <- data
LUAD_clin_snv <- clin_snv

LUSC_data <- data
LUSC_clin_snv <-clin_snv

# 提取 RiskGroup 为 High_Risk 的样本
high_risk_samples <- CoxData %>% filter(RiskGroup == "High_Risk") %>% pull(Sample)

low_risk_samples <- CoxData %>% filter(RiskGroup == "Low_Risk") %>% pull(Sample)

LUAD_data <-subset(LUAD_data,LUAD_data$Tumor_Sample_Barcode %in% low_risk_samples)

LUAD_maf <- read.maf(LUAD_datah)

pdf('./result/figure/12_Fig12C_4_LUSC_high.pdf',width = 9,height = 12)
maftools::oncoplot(LUAD_maf)
dev.off()

　class(LUAD_maf)
pdf('./result/figure/12_Fig12D_2_LUAD_low.pdf',width = 8,height = 10)

plotmafSummary(maf = LUAD_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

dev.off()


Prognostic_gene<-readRDS('./result/Prognostic_gene.rds')
cat(Prognostic_gene, sep = " ")
