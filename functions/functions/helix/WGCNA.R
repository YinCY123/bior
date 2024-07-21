
####WGCNA####
library(WGCNA)
library(reshape2)
library(stringr)

dir.create('WGCNA')
options(stringsAsFactors = FALSE)

type <- "unsigned"
corType <- "pearson"
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

dataExpr <- exp  
dim(dataExpr)

if(continuous){
  traitData=as.matrix(pheno_score)
}else {
  traitData=model.matrix(~0+ group)
  colnames(traitData)=gsub("group","",colnames(traitData))
  rownames(traitData)=names(group)
}



#数据筛选#
#选择挑选方差前5000的基因进行后续分析
m.mad <- apply(dataExpr,1,var)    
if(nrow(dataExpr)>5000){
  dataExprVar <- dataExpr[order(m.mad,decreasing = T)[1:5000],]
}else{
  dataExprVar <- dataExpr[which(m.mad >
                                  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
}

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))   #方差前5000的基因

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK          #TRUE
#如果gsg$allOK的结果为TRUE，证明没有缺失值，可以直接下一步。如果为FALSE，则需要用以下函数进行删除缺失值。

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


dir.create("WGCNA/plot")
dir.create("WGCNA/data")
dir.create("WGCNA/Table")
WGCNA_list=list()
#软阈值筛选#
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)
pdf("WGCNA/plot/软阈值筛选.pdf",height = 8/2.54,width = 12/2.54)
par(mfrow = c(1,2),mgp=c(1, 0.2, 0),mar=c(2,2,1,0.5),ps=6,tck = -0.02)


# 横轴是Soft threshold(power),纵轴是无标度网络的评估参数,数值越高,网络越符合无标度特征(non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
abline(h=0.85,col="red")# 筛选标准 R-square=0.85

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")

dev.off()
power = sft$powerEstimate   #设置软阈值beta值

WGCNA_list[["sft_fitIndices"]]=sft$fitIndices

##无满足条件的power时选用经验power
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

WGCNA_list[["power"]]=power
#网络构建#
net = blockwiseModules(dataExpr, power = power, 
                       maxBlockSize = nGenes,#计算机能处理的最大模块的基因数量
                       TOMType = type, minModuleSize = 30,     #模块数量太多时，提高模块内最低基因数值
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25, #合并模块的阈值，越大模块越少
                       numericLabels = TRUE, #返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
                       pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("exprMat.tom"),
                       verbose = 3)
table(net$colors)#0(grey)表示未分入任何模块的基因      

f=0
while (ncol(net$MEs)>16){
  if(f<2){f=f+1}else{f=f+5}
  net = blockwiseModules(dataExpr, power = power, 
                         maxBlockSize = nGenes,#计算机能处理的最大模块的基因数量
                         TOMType = type, minModuleSize = 30+10*f,     #模块数量太多时，提高模块内最低基因数值
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25, #合并模块的阈值，越大模块越少
                         numericLabels = TRUE, #返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
                         pamRespectsDendro = FALSE,
                         saveTOMs=TRUE, corType = corType,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = paste0("exprMat.tom"),
                         verbose = 3)
}


WGCNA_list[["minModuleSize"]]=30+10*f
WGCNA_list[["net"]]=net


#模块可视化#
##层级聚类树
moduleLabels = net$colors     
moduleColors = labels2colors(moduleLabels)
WGCNA_list[["moduleColors"]]=moduleColors
WGCNA_list[["moduleLabels_count"]]=data.frame(table(moduleColors))
pdf("WGCNA/plot/基因聚类模块.pdf",height = 6,width = 10)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()
##绘制模块之间相关性图
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf("WGCNA/plot/模块间的相关性热图.pdf",height = 6,width = 6)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(1,3,2,4),
                      marHeatmap = c(3,3,1,2), # marDendro/marHeatmap设置下左上右的边距
                      plotDendrograms = T,
                      xLabelsAngle = 90)
dev.off()

# #可视化基因网络 (TOM plot)#
# load(net$TOMFiles[1], verbose=T)
# TOM <- as.matrix(TOM)
# dissTOM = 1-TOM
# #抽样
# nSelect <- 400
# select <- sample(nGenes, size = nSelect)
# selectTOM <- dissTOM[select, select]
# 
# # 根据TOM矩阵重新对基因聚类；
# selectTree <- hclust(as.dist(selectTOM), method = "average")
# #提取选择基因的模块颜色；
# selectColors <- moduleColors[select]
# # 对TOM矩阵进行指数转化，使热图能展示更丰富的信息；
# plotDiss <- selectTOM^7
# diag(plotDiss) <- NA
# #绘制热图，默认数值越大（越接近1）颜色越深（底色）；
# pdf("WGCNA/plot/基因网络TOM.pdf",height = 6,width = 6)
# TOMplot((1-plotDiss), selectTree, selectColors,
#         main = "Network heatmap plot, selected genes")
# dev.off()

### 模块与表型数据关联
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## 模块与表型相关性热图
textMatrix = paste(signif(modTraitCor, 4), "\n(", sprintf("%.0e",signif(modTraitP, 1)), ")", 
                   sep = "")# signif表示保留几位小数
dim(textMatrix) = dim(modTraitCor)
pdf("WGCNA/plot/模块和性状的关系.pdf",width = 8,height = 8)
par(mar = c(3,6,4,3))
labeledHeatmap(Matrix = modTraitCor, 
               xLabels = colnames(traitData),xLabelsAngle = 0,
               yLabels = colnames(MEs_col),
               cex.lab = 0.8,
               ySymbols = colnames(MEs_col), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.8, textAdj = c(0.5, 0.5),
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

#把所有参数保存到一个list中
WGCNA_list[["labeledHeatmap_data"]]=list(
  Matrix = modTraitCor,
  xLabels = colnames(traitData),
  xLabelsAngle = 0,
  yLabels = colnames(MEs_col),
  cex.lab = 0.8,
  ySymbols = colnames(MEs_col), colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix, setStdMargins = FALSE,
  cex.text = 0.8, textAdj = c(0.5, 0.5),
  zlim = c(-1,1),
  main = paste("Module-trait relationships")
)
#另外写一个用刚刚保存的参数进行绘图的函数
WGCNA_list[["labeledHeatmap"]]=function(xLabelsAngle = WGCNA_list$labeledHeatmap_data$xLabelsAngle){
  WGCNA::labeledHeatmap(Matrix = WGCNA_list$labeledHeatmap_data$Matrix, 
                 xLabels = WGCNA_list$labeledHeatmap_data$xLabels,xLabelsAngle = xLabelsAngle,
                 yLabels = WGCNA_list$labeledHeatmap_data$yLabels,
                 cex.lab = WGCNA_list$labeledHeatmap_data$cex.lab,
                 ySymbols = WGCNA_list$labeledHeatmap_data$ySymbols, colorLabels = WGCNA_list$labeledHeatmap_data$colorLabels,
                 colors = WGCNA_list$labeledHeatmap_data$colors,
                 textMatrix = WGCNA_list$labeledHeatmap_data$textMatrix, setStdMargins = WGCNA_list$labeledHeatmap_data$setStdMargins,
                 cex.text = WGCNA_list$labeledHeatmap_data$cex.text, textAdj = WGCNA_list$labeledHeatmap_data$textAdj,
                 zlim = WGCNA_list$labeledHeatmap_data$zlim,
                 main = WGCNA_list$labeledHeatmap_data$main)
}


## 计算模块与基因的相关性矩阵
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}
WGCNA_list[["geneModuleMembership"]]=geneModuleMembership

## 计算性状与基因的相关性矩阵
if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}
WGCNA_list[["geneTraitCor"]]=geneTraitCor

## 两个相关性矩阵联合起来,指定感兴趣模块进行分析
WGCNA_list[["modTraitCor"]]=modTraitCor
WGCNA_list[["modTraitP"]]=modTraitP

write.csv(modTraitCor,'WGCNA/Table/modTraitCor.csv')
write.csv(modTraitP,'WGCNA/Table/modTraitP.csv')

#筛选关联性最强的模块名字
#先找出p值大于0.05的模块,再找出关联值最大的模块名字
module_names=data.frame()
if(!(exists("case")|is.null(case))){
  for (i in 1:nrow(modTraitCor)){
    if (modTraitP[i,case]<0.05){
      module_names=rbind(module_names,
                         data.frame(row.names=rownames(modTraitCor)[i],modCor=signif(modTraitCor[i,case], 4),case=case))
    }
  }
}else{
  for (i in 1:nrow(modTraitCor)){
    temp=modTraitCor[i,]
    j=which(abs(temp)==max(max(abs(temp))))[1]
    if (modTraitP[i,j]<0.05){
      module_names=rbind(module_names,data.frame(row.names=rownames(modTraitCor)[i],modCor=signif(modTraitCor[i,j], 4),case=colnames(modTraitCor)[j]))
    }
  }
}

module_top=module_names%>%slice_max(n = 2,order_by = abs(modCor),with_ties = F)
WGCNA_list[["module_top"]]=module_top
WGCNA_list[["module_names"]]=module_names

####相关性绝对值最大的模块####
module=rownames(module_top)[1]
module=gsub('ME','',module)         #module_trait relationship 中感兴趣的模块名字
pheno = module_top[1,2]            #module_trait relationship 中感兴趣的性状名字
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

pdf("WGCNA/plot/模块中基因相关性散点图1.pdf",height = 6,width = 6)
#sizeGrWindow(7, 7)
par(mfrow = c(1,1),mar=c(6,6,6,6))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#提取指定模块内的基因名#
gene <- rownames(geneTraitCor)[moduleGenes]    
write.csv(gene,"WGCNA/Table/Module_genes_top1.csv")

if(nrow(module_top)==2){
  ####相关性绝对值top2的模块####
  module=rownames(module_top)[2]
  module=gsub('ME','',module)         #module_trait relationship 中感兴趣的模块名字
  pheno = module_top[2,2]            #module_trait relationship 中感兴趣的性状名字
  modNames = substring(colnames(MEs_col), 3)
  # 获取关注的列
  module_column = match(module, modNames)
  pheno_column = match(pheno,colnames(traitData))
  # 获取模块内的基因
  moduleGenes = moduleColors == module
  
  pdf("WGCNA/plot/模块中基因相关性散点图2.pdf",height = 6,width = 6)
  #sizeGrWindow(7, 7)
  par(mfrow = c(1,1),mar=c(6,6,6,6))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                     abs(geneTraitCor[moduleGenes, pheno_column]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for", pheno),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  #提取指定模块内的基因名#
  gene2 <- rownames(geneTraitCor)[moduleGenes]    
  write.csv(c(gene,gene2),"WGCNA/Table/Module_genes_top2.csv")
}

####提取自选模块基因####
WGCNA_list[["moduel_select"]]<-function(module){
  pheno=WGCNA_list$module_names[paste0('ME',module),2]
  modNames1 = substring(rownames(WGCNA_list$modTraitCor), 3)
  module_column = match(module, modNames)
  pheno_column = match(pheno,colnames(WGCNA_list$modTraitCor))
  moduleGenes = WGCNA_list$moduleColors == module
  
  pdf("WGCNA/plot/模块中基因相关性散点图use.pdf",height = 6,width = 6)
  #sizeGrWindow(7, 7)
  par(mfrow = c(1,1),mar=c(6,6,6,6))
  verboseScatterplot(abs(WGCNA_list$geneModuleMembership[moduleGenes, module_column]),
                     abs(WGCNA_list$geneTraitCor[moduleGenes, pheno_column]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for", pheno),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  #提取指定模块内的基因名#
  gene_use <- rownames(WGCNA_list$geneTraitCor)[moduleGenes]    
  write.csv(gene_use,"WGCNA/Table/Module_gene_use.csv")
}

save(WGCNA_list,file = "WGCNA/data/WGCNA_list.rda")
