setwd("/data/home/jiangminghui/NSCLC")
library(WGCNA)
library(reshape2)
library(stringr)
library(Seurat)
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
library(DESeq2)
library(edgeR)
library(igraph)
library(tinyarray)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(GSVA)
library(paletteer)
library(survival)
library(survminer)
library(cluster)
library(sva)

dim(expr)
dim(clinical_data)
dim(TCGA_NSCLC_clin)
############################################################
#
#-------------------------读取数据---------------------
#
############################################################
TCGA_NSCLC_clin<-readRDS('./data/02_TCGA_NSCLC_clin.rds')
expr<- readRDS('./data/02_TCGA_expr.rds')
dim(expr)

type = "unsigned"
corType = "pearson"

corFnc = if (corType == "pearson") {
  corFnc = cor
} else {
  corFnc = bicor
}

maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)

range(expr)
expr1 <- log2(expr+1)
range(expr1)

dim(expr1)

m.mad <- apply(expr1,1,mad,na.rm=T)

dataExprVar <- expr1[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]



## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)
table(gsg$goodSamples)
gsg$allOK


if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)
head(dataExpr)[,1:8]
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
abline(h = 180, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = 180, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
dataExpr = dataExpr[keepSamples, ]
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
pdf('./result/figure/03_WGCNA_Soft_threshold_power.pdf',width = 8,height = 6)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.87,col="red")
dev.off()


power = sft$powerEstimate
power
pdf('./result/figure/03_WGCNA_mean_Connectivity.pdf',width = 8,height = 6 )
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
dev.off()

power = sft$powerEstimate
power



net_0807= blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0('NSCLC', ".tom"),
                       verbose = 3)
saveRDS(net_0807,file = './result/data/03_WGCNA_net.rds')

#net <-readRDS('./result/data/03_WGCNA_net.rds')
net <-net_0807


table(net_0807$colors)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
pdf('./result/figure/03_Fig3E_WGCNA_Cluster_tree.pdf',width = 8,height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()


##绘制模块之间的相关性图
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
pdf('./result/figure/03_Sup_Fig_WGCNA_Correlation_modules.pdf',height = 10,width = 8)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

gsva_mat <- as.data.frame(t(gsva(as.matrix(expr1),gset.idx.list = list(MCU_Geneset),method = 'ssgsea'))) #计算得分
gsva_mat$Group <- ifelse(gsva_mat$V1 > median(gsva_mat$V1), "High", "Low")

head(common_samples)

common_samples <- intersect(rownames(MEs_col), rownames(gsva_mat))


MEs_col_filtered <- MEs_col[common_samples, ]
gsva_mat_filtered <- gsva_mat[common_samples, ]

dim(MEs_col_filtered)

dim(gsva_mat_filtered)

MEs_colpheno <- orderMEs(cbind(MEs_col_filtered, gsva_mat_filtered))

saveRDS(gsva_mat_filtered,file = './result/data/03_WGCNA_gsva_mat_filtered.rds')

load(net$TOMFiles[1], verbose=T)

#TOM <- as.matrix(TOM)

#dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
#plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
#diag(plotTOM) = NA
# Call the plot function

# 这一部分特别耗时，行列同时做层级聚类

#pdf('./result/figure/03_Fig3E_WGCNA_TOM_plot.pdf',width = 8,height = 6)
#TOMplot(plotTOM, net$dendrograms, moduleColors, 
#        main = "Network heatmap plot, all genes")
#dev.off()


str(dataExpr)

#明确样本数和基因数
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
#首先针对样本做个系统聚类树
datExpr_tree<-hclust(dist(dataExpr), method = "average")

saveRDS(datExpr_tree,file = './result/data/03_WGCNA_datExpr_tree.rds')

par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)

str(sample_colors)

head(gsva_mat_filtered)
sample_colors <- numbers2colors(as.numeric(factor(gsva_mat_filtered$Group)), 
                                colors = c("white","blue"),signed = FALSE)

pdf('./result/figure/03_Fig3B_Sample_Cluster_Map.pdf',width = 20,height = 8)
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

dev.off()

head(gsva_mat_filtered)
str(gsva_mat)
design=model.matrix(~0+ gsva_mat_filtered$Group)
colnames(design)=levels(gsva_mat_filtered$Group)
moduleColors <- labels2colors(net$colors)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes

saveRDS(MEs0,file = './result/data/03_WGCNA_MEs0.rds')

MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


head(moduleTraitCor)

str(design)
str(gsva_mat_filtered)

par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
pdf('./result/figure/03_Sup_Fig_WGCNA_Module-trait_relationships.pdf',width = 10,height = 20)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c('group','score'),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

str(moduleTraitCor)




sampleName = rownames(dataExpr)
dim(dataExpr) #1030 31846
traitData = gsva_mat[match(sampleName, rownames(gsva_mat)), ]
# Convert Group to a factor
traitData$Group <- as.factor(traitData$Group)

# Relevel the factor so that "High" is 1 and "Low" is 2
traitData$Group <- relevel(traitData$Group, ref = "Low")

# Convert the factor to numeric
traitData$Group <- as.numeric(traitData$Group)
colnames(traitData) <- c('MCU-RGs Score','Group')
table(traitData$Group)
#traitData$Group <-as.numeric(as.factor(traitData$Group))
dim(traitData) #1030    2
#traitData$HRD_class <- ifelse(traitData$HRD>46,1,0)
#traitData <- as.numeric(as.matrix(traitData))
head(traitData)
### 模块与表型数据关联
str(MEs_col) #1030
if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}



pdf('./result/figure/03_Sup_Fig_WGCNA_Module-Score_heatmap.pdf',width = 8,height =8 )
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.3, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.3, zlim = c(-1,1),
               main = paste("Module-Score Heatmap"))
dev.off()


# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
#module = c('magenta','blue','green')
module = c('blue')
#module = "grey60"
pheno = "MCU-RGs Score"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

table(moduleGenes)
if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf('./result/figure/03_Fig3G_WGCNA_Scatterplot_GS_MM_01.pdf',width = 6,height = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h=0.4,v = 0.6,col="red")


module_Genes <- rownames(geneModuleMembership[moduleGenes,])

print(module)
print(pheno)

dim(geneModuleMembership[moduleGenes, ])
dim(geneTraitCor[moduleGenes, ])

module_Genes <- rownames(geneModuleMembership[moduleGenes,])
length(module_Genes)


MM <- abs(geneModuleMembership[moduleGenes, module_column])
GS <- abs(geneTraitCor[moduleGenes, pheno_column])

save(MM,GS,file = './result/data/03_WGCNA_MM_GS.RData')

WGCNA_keyGene <- names(MM)[which(MM>0.6&GS>0.4)]
length(WGCNA_keyGene)
save(WGCNA_keyGene,file = './result/data/WGCNA_keyGene.RData')




# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
#module = c('magenta','blue','green')
module = c('blue')
#module = "grey60"
pheno = "MCU-RGs Score"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

table(moduleGenes)
if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf('./result/figure/03_Fig3G_WGCNA_Scatterplot_GS_MM_02.pdf',width = 6,height = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h=0.75,v = 0.8,col="red")


module_Genes <- rownames(geneModuleMembership[moduleGenes,])

print(module)
print(pheno)

dim(geneModuleMembership[moduleGenes, module_column])
dim(geneTraitCor[moduleGenes, pheno_column])

module_Genes <- rownames(geneModuleMembership[moduleGenes,])
length(module_Genes)

save(MM,GS,file = './result/data/03_WGCNA_MM_GS_02.RData')

MM <- abs(geneModuleMembership[moduleGenes, module_column])
GS <- abs(geneTraitCor[moduleGenes, pheno_column])

WGCNA_keyGene2 <- names(MM)[which(MM>0.6&GS>0.2)]
length(WGCNA_keyGene2)
save(WGCNA_keyGene2,file = './result/data/WGCNA_keyGene2.RData')


# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
#module = c('magenta','blue','green')
module = c('green')
#module = "grey60"
pheno = "MCU-RGs Score"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

table(moduleGenes)
if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf('./result/figure/03_Fig3G_WGCNA_Scatterplot_GS_MM_03.pdf',width = 6,height = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h=0.75,v = 0.8,col="red")


module_Genes <- rownames(geneModuleMembership[moduleGenes,])

print(module)
print(pheno)

dim(geneModuleMembership[moduleGenes, module_column])
dim(geneTraitCor[moduleGenes, pheno_column])

module_Genes <- rownames(geneModuleMembership[moduleGenes,])
length(module_Genes)

save(MM,GS,file = './result/data/03_WGCNA_MM_GS_03.RData')

MM <- abs(geneModuleMembership[moduleGenes, module_column])
GS <- abs(geneTraitCor[moduleGenes, pheno_column])

WGCNA_keyGene3 <- names(MM)[which(MM>0.6&GS>0.2)]
length(WGCNA_keyGene3)
save(WGCNA_keyGene3,file = './result/data/WGCNA_keyGene3.RData')


















