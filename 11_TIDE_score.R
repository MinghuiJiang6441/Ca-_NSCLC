library("ggstatsplot")
library(magrittr)

TIDE_Res <-read.csv('./data/TIDE_Res.csv')

head(TIDE_Res,3)

TIDE_Res$Sample <- colnames(expr)[TIDE_Res$Patient]

rownames(TIDE_Res) <- TIDE_Res$Sample 
TIDE_Res <-TIDE_Res[CoxData$Sample,]
TIDE_Res$group <- CoxData$RiskGroup

saveRDS(TIDE_Res,'./result/data/TIDE_Res.rds')
TIDE_Res<-readRDS('./result/data/11_TIDE_Res.rds')

p1 <- ggbetweenstats(
  data  = TIDE_Res,
  x     = group,
  y     = Exclusion,
  plot.type = "box",#图标类型，不加参数则默认为boxviolin，其它可选为：box、violin
  title = "Distribution of Exclusion across Risk group"
)

p1#上述为ggbetweenstats()的最简单函数调用

p2 <- ggbetweenstats(
  data  = TIDE_Res,
  x     = group,
  y     = Dysfunction,
  plot.type = "box",#图标类型，不加参数则默认为boxviolin，其它可选为：box、violin
  title = "Distribution of Dysfunction across Risk group"
)

p2#


p3 <- ggbetweenstats(
  data  = TIDE_Res,
  x     = group,
  y     = TIDE,
  plot.type = "box",#图标类型，不加参数则默认为boxviolin，其它可选为：box、violin
  title = "Distribution of TIDE across Risk group"
)

p3#
pdf('./result/figure/11_Fig11E_TIDE_score.pdf',width = 12,height = 4,onefile = F)
p1 + p2 +p3

dev.off()

DATA
fit <- survfit(Surv(A7_OS.time , A7_OS.event) ~ RiskGroup, data = TIDE_Res)

TIDE_Res$TIDE_Group <- ifelse(TIDE_Res$TIDE > median(TIDE_Res$TIDE), "HTIDE", "LTIDE")

TIDE_Res$Risk_TIDE_Group <- paste(TIDE_Res$group, TIDE_Res$TIDE_Group, sep = "-")


table(TIDE_Res$Risk_TIDE_Group )


CoxData$Risk_TIDE_Group <- TIDE_Res$Risk_TIDE_Group

fit <- survfit(Surv(OS.time , OS) ~ Risk_TIDE_Group, data = CoxData)

km<-ggsurvplot(
  fit, 
  data = CoxData, 
  size = 1,                 # change line size
  palette = 
    c("#E7B800", "#2E9FDF",'#FFCB9A','#0F6466'),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
pdf('./result/figure/11_Fig11F_km_curve.pdf',width = 6,height = 6,onefile = F)

km
dev.off()

############################################################
#
#--------------TMB————————————————————————————————
#
############################################################

LUAD_maf <- read.maf(maf = "data/Somatic_LUAD.MAF")

LUAD_tmb_table <- tmb(maf = LUAD_maf, logScale = FALSE)

LUSC_maf <- read.maf(maf = "data/LUSC_ Somatic.MAF")

LUSC_tmb_table <- tmb(maf = LUSC_maf, logScale = FALSE)

tmb_table <-rbind(LUAD_tmb_table,LUSC_tmb_table) %>%as.data.frame()

head(CoxData,3)


rownames(tmb_table) <-substr(tmb_table$Tumor_Sample_Barcode,start = 1,stop = 12)

tmb_table<-tmb_table[CoxData$Sample,]

tmb_table$total_perMB[is.na(tmb_table$total_perMB)] <- 0

tmb_table$TMB_group <- ifelse(tmb_table$total_perMB > median(tmb_table$total_perMB), "HTMB", "LTMB")

CoxData$TMB_group <-tmb_table$TMB_group

CoxData$Risk_TMB_Group <- paste(CoxData$RiskGroup, CoxData$TMB_group, sep = "-")

saveRDS(object = tmb_table,'./data/tmb_table.rds')
saveRDS(CoxData,file = './result/data/CoxData.rds')
high_risk_prop <- prop.table(table(CoxData$RiskGroup))
low_risk_prop <- prop.table(table(low_risk_tmb$group))


risk_group <- CoxData$RiskGroup
tmb_group <- CoxData$TMB_group

# 计算每组的相对比例
prop_table <- prop.table(table(risk_group, tmb_group), margin = 1)
plot_data <- as.data.frame(prop_table)
colnames(plot_data) <- c("RiskGroup", "TMBGroup", "Proportion")


p <-ggplot(plot_data, aes(x = RiskGroup, y = Proportion, fill = TMBGroup)) +
  scale_fill_manual(values = c(  "#D91A2A","#03A6A6" )) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Relative Proportion of High and Low TMB Scores in High and Low Risk Groups",
       x = "Risk Group",
       y = "Proportion") +
  theme_minimal()
ggsave(filename = './result/figure/11_Fig11G_barplot.pdf',p,width = 5,height = 6)



head(CoxData,2)


fit <- survfit(Surv(OS.time , OS) ~ Risk_TMB_Group, data = CoxData)

km<-ggsurvplot(
  fit, 
  data = CoxData, 
  size = 1,                 # change line size
  palette = 
    c("#E7B800", "#2E9FDF",'#FFCB9A','#0F6466'),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

pdf('./result/figure/11_Fig11H_KM_curve.pdf',width = 6,height = 6,onefile = F)
km
dev.off()