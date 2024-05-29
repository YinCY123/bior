####WGCNA####
dir.create('WGCNA')
exp1=as.matrix(combat_expr)

options(stringsAsFactors = FALSE)

type <- "unsigned"
corType <- "pearson"
corFnc = ifelse(corType=="pearson", WGCNA::cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

#导入表型数据#
#ssGSEA分组打分
library(GSVA)
min.sz = 2 # 单个基因集基因数下限
max.sz = 10000 # 基因数上限
parallel.sz= 50 # 并行处理cpu数量
mx.diff= T # ES 计算为最大正和负随机游走偏差之间的幅度差。
tau=0.25 # 默认1，tau=1 when method="gsva" and tau=0.25 when method="ssgsea" 
method='ssgsea' # 默认
kcdf = "Gaussian" # 默认"Gaussian"，当输入数值为整数时，设置为"Poisson"

ssgsea = gsva(expr=exp1, gset.idx.list=pheno_symbol, 
              method=method, kcdf=kcdf, min.sz=min.sz, max.sz=max.sz, 
              parallel.sz=parallel.sz, mx.diff=mx.diff,tau=0.25 )

write.csv(ssgsea,'WGCNA/ssGSEA.csv')
ssgsea_group=t(ssgsea)
colnames(ssgsea_group)=phenotype

dataExpr <- exp1  
dim(dataExpr)        
datTraits <- data.frame(gsm=rownames(ssgsea_group),group=ssgsea_group) #group是表型打分用作分组信息
colnames(datTraits)=c('gsm','group')

#数据筛选#
#选择挑选方差前5000的基因进行后续分析
m.mad <- apply(dataExpr,1,var)    
dataExprVar <- dataExpr[order(m.mad,decreasing = T)[1:5000],]

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
head(dataExpr)[,1:8]

## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                colors = c("blue","red"),signed = FALSE)
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(sampleTree, sample_colors,
                    groupLabels = "Group",
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

#软阈值筛选#
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)
pdf("WGCNA/软阈值筛选.pdf",height = 5.5,width = 8)
par(mfrow = c(1,2))
cex1 = 0.9

# 横轴是Soft threshold(power),纵轴是无标度网络的评估参数,数值越高,网络越符合无标度特征(non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")# 筛选标准 R-square=0.85

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")

while (!is.null(dev.list()))  dev.off()
power = sft$powerEstimate   #设置软阈值beta值
power

write.csv(sft$fitIndices,"WGCNA/sft_fitIndices.csv",row.names = F,quote = F)

##无满足条件的power时选用经验power
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

#网络构建#
cor=WGCNA::cor
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
  f=f+1
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


write.csv(data.frame(table(net$colors)),'WGCNA/Module_Nofgenes.csv')
saveRDS(net,"WGCNA/net.rds")

cor=stats::cor

#模块可视化#
##层级聚类树
pdf("WGCNA/基因聚类模块.pdf",height = 6,width = 10)
moduleLabels = net$colors     
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,cex.axis=0.6,cex.lab=1,cex.main=1)

while (!is.null(dev.list()))  dev.off()
##绘制模块之间相关性图
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf("WGCNA/模块间的相关性热图.pdf",height = 6,width = 6)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marHeatmap = c(3,3,1,2), # marDendro/marHeatmap设置下左上右的边距
                      plotDendrograms = F,
                      xLabelsAngle = 90,cex.axis=0.6,cex.lab=1,cex.main=1)
while (!is.null(dev.list()))  dev.off()

#可视化基因网络 (TOM plot)#
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
#抽样
nSelect <- 400
set.seed(10)
select <- sample(nGenes, size = nSelect)
head(select)
selectTOM <- dissTOM[select, select]
selectTOM[1:6,1:6]

# 根据TOM矩阵重新对基因聚类；
selectTree <- hclust(as.dist(selectTOM), method = "average")
#提取选择基因的模块颜色；
selectColors <- moduleColors[select]
# 对TOM矩阵进行指数转化，使热图能展示更丰富的信息；
plotDiss <- selectTOM^7
diag(plotDiss) <- NA
#绘制热图，默认数值越大（越接近1）颜色越深（底色）；
pdf("WGCNA/基因网络TOM.pdf",height = 6,width = 6)
TOMplot((1-plotDiss), selectTree, selectColors,
        main = "Network heatmap plot, selected genes")
while (!is.null(dev.list()))  dev.off()

#关联表型数据#
traitData=model.matrix(~0+ datTraits$group)
colnames(traitData)=phenotype

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
textMatrix = paste(signif(modTraitCor, 4), "\n(", signif(modTraitP, 1), ")", 
                   sep = "")# signif表示保留几位小数
dim(textMatrix) = dim(modTraitCor)
pdf("WGCNA/模块和性状的关系.pdf",width = 8,height = 8)
par(mar = c(3,6,4,3))
labeledHeatmap(Matrix = modTraitCor, 
               xLabels = colnames(traitData),xLabelsAngle = 0,
               yLabels = colnames(MEs_col),
               cex.lab = 0.6,cex.main=1,
               ySymbols = colnames(MEs_col), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.2, textAdj = c(0.5, 0.5),
               zlim = c(-1,1),mar = c(0,0,0,0),
               main = paste("Module-trait relationships"))

while (!is.null(dev.list()))  dev.off()


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

## 两个相关性矩阵联合起来,指定感兴趣模块进行分析
modTraitCor
modTraitP
write.csv(modTraitCor,'WGCNA/modTraitCor.csv')
write.csv(modTraitP,'WGCNA/modTraitP.csv')

#筛选关联性最强的模块名字
#先找出p值大于0.05的模块,再找出关联值最大的模块名字
module_names=data.frame()
for (c in 1:nrow(modTraitCor)){
  if (modTraitP[c]<0.05){
    module_names=rbind(module_names,data.frame(row.names=rownames(modTraitCor)[c],modCor=modTraitCor[c]))
  }
}
module=rownames(module_names)[which(abs(module_names$modCor)==max(abs(module_names$modCor)))]
module=gsub('ME','',module)         #module_trait relationship 中感兴趣的模块名字
pheno = phenotype            #module_trait relationship 中感兴趣的性状名字
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

pdf("WGCNA/模块中基因相关性散点图.pdf",height = 6,width = 6)
#sizeGrWindow(7, 7)
par(mfrow = c(1,1),mar=c(6,6,6,6))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
while (!is.null(dev.list()))  dev.off()

#提取指定模块内的基因名#
gene <- rownames(geneTraitCor)[moduleGenes]    
write.table(gene,"WGCNA/Module_genes.txt",sep='\t',quote = F,row.names = F,col.names = F)
write.csv(gene,"WGCNA/Module_genes.csv")

#方差前5000筛选的与表型最相关模块内的基因，与差异基因取交集
gene1=dplyr::intersect(gene,deg_gene)

write.table(gene1,paste("WGCNA/","WGCNA",disease,phenotype,"gene1.txt",collapse = "_"),sep='\t',quote = F,row.names = F,col.names = F)
write.csv(gene1,paste("WGCNA/","WGCNA",disease,phenotype,"gene1.csv",collapse = "_"))