
####参数检测####
if(class(exp)[1]!="matrix")exp=as.matrix(exp)
if(length(unique(group))!=2)warning('分组不是2组，请检查')
if(!(case%in%group&control%in%group))stop('case和control与group中的信息未对应上')
if(!all(colnames(exp)==names(group)))stop('样本名未对应上')
if(sum(is.na(exp))>=1)stop('表达矩阵有空值')
if(is.null(colors))stop("缺少颜色参数")
####差异分析####
dir.create("DEG")
setwd("DEG/")
if(is_count){
  cat("你选择了使用count矩阵计算差异基因，这会非常耗时，另外请检查数据是否为count矩阵")
if(max(table(group))>15){  
  ####edgeR####
  library(edgeR)
  
  group_deg <- factor(group,levels = c(control,case))
  design <- model.matrix(~0+group_deg)
  row.names(design) <- colnames(exp)
  
  exp[exp<0] <- 0#负数设为0
  DGElist <- DGEList(counts = exp,group = group_deg)
  ## 过滤基因
  keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2#至少有两个样本cpm大于1的基因
  
  DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
  ## 标准化
  DGElist <- calcNormFactors( DGElist )
  DGElist <- estimateGLMCommonDisp(DGElist, design)
  DGElist <- estimateGLMTrendedDisp(DGElist, design)
  DGElist <- estimateGLMTagwiseDisp(DGElist, design)
  
  fit <- glmFit(DGElist, design)
  results <- glmLRT(fit, contrast = c(-1, 1)) 
  nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
  nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
  head(nrDEG_edgeR)
  
  dif <- na.omit(nrDEG_edgeR)
  dif$gene_symbol <- row.names(dif)
  library(dplyr)
  if(P=='p'){
    dif.up <- dif %>%
      dplyr::filter(logFC > logFC_value & PValue < P_value)
    dif.down <- dif %>%
      dplyr::filter(logFC < (-logFC_value) & PValue < P_value)
    dif <- dif %>% dplyr::mutate(group = dplyr::case_when(
      gene_symbol %in% dif.up$gene_symbol ~ 'up',
      gene_symbol %in% dif.down$gene_symbol ~ 'down',
      TRUE ~ "no"
    ))}else{
      dif.up <- dif %>%
        dplyr::filter(logFC > logFC_value & FDR < P_value)
      dif.down <- dif %>%
        dplyr::filter(logFC < (-logFC_value) & FDR < P_value)
      dif <- dif %>% dplyr::mutate(group = dplyr::case_when(
        gene_symbol %in% dif.up$gene_symbol ~ 'up',
        gene_symbol %in% dif.down$gene_symbol ~ 'down',
        TRUE ~ "no"
      ))
    }
  
  DEG=data.frame(name=dif$gene_symbol,
                 logFC=dif$logFC,
                 P.value=dif$PValue,
                 P.adj=dif$FDR,
                 group=dif$group)
  ####输出数据####
  dir.create("data")
  dir.create("Table")
  print(table(DEG$group))
  
  DEG_list<-list(method="edgeR",
                 case=case,
                 control=control,
                 DEG=DEG,P=P,
                 logFC_value=logFC_value,P_value=P_value,
                 feature_type=feature_type,colors=colors)
  save(DEG_list,file = "data/DEG_list.rda")
  write.csv(DEG,'Table/DEG.csv')
  
  exp_plot=log2(cpm(DGElist)+1)#用于后续绘图
}else{
  ####DEseq2####
  library(DESeq2)
  group_deg <- factor(group,levels = c(control,case))
  # 创建一个数据框用于存储样本的分组信息，行名为样本名，列名为分组信息
  colData <- as.data.frame(group_deg)
  # DESeq2 要求输入数据是由整数组成的矩阵，且没有经过标准化
  exp_int <- exp
  exp_int <- apply(exp_int, 2, as.integer)
  rownames(exp_int) <- rownames(exp)
  
  # 构建DESeqDataSet对象，也就是dds矩阵，将基因计数数据、样本分组信息和设计矩阵关联起来
  dds <- DESeqDataSetFromMatrix(countData = exp_int, # 表达矩阵
                                colData = colData,        # 表达矩阵列名和分组信息的对应关系
                                design = ~ group_deg)         # group_deg为colData中的group_deg，也就是分组信息
  
  # 进行差异表达分析
  dds <- DESeq(dds)
  
  # 提取差异表达结果，进行对比，这里contrast参数指定了对比的组别
  # contrast参数必须写成下面三个元素的向量格式，且顺序不能反
  res <- results(dds, contrast = c("group_deg", case,control))
  
  # 按照padj（调整后的p值）的大小对差异结果进行排序（只有DESeq2需要，limma和edgeR会自动排好）
  resOrdered <- res[order(res$padj), ]
  
  # 将差异表达结果转换为数据框
  resOrdered <- as.data.frame(resOrdered)
  dif <- na.omit(resOrdered)
  dif$gene_symbol <- row.names(dif)
  if(P=="p.adj"){
    dif.up <- dif %>%
      dplyr::filter(log2FoldChange > logFC_value & padj < P_value)#筛选上调基因，logFC可调整，P值不能大于0.05
    dif.down <- dif %>%
      dplyr::filter(log2FoldChange < (-logFC_value) & padj < P_value)#筛选下调基因，logFC可调整，P值不能大于0.05
  }else{
    dif.up <- dif %>%
      dplyr::filter(log2FoldChange > logFC_value & pvalue < P_value)#筛选上调基因，logFC可调整，P值不能大于0.05
    dif.down <- dif %>%
      dplyr::filter(log2FoldChange < (-logFC_value) & pvalue < P_value)#筛选下调基因，logFC可调整，P值不能大于0.05
  }
  #添加差异基因信息
  dif <- dif %>% dplyr::mutate(group = dplyr::case_when(
    gene_symbol %in% dif.up$gene_symbol ~ "up",
    gene_symbol %in% dif.down$gene_symbol ~ "down",
    TRUE ~ "no"
  ))
  
  
  DEG=data.frame(name=dif$gene_symbol,
                 logFC=dif$log2FoldChange,
                 P.value=dif$pvalue,
                 P.adj=dif$padj,
                 group=dif$group)
  
  ####输出数据####
  dir.create("data")
  dir.create("Table")
  print(table(DEG$group))
  
  DEG_list<-list(method="DESeq2",
                 case=case,
                 control=control,
                 DEG=DEG,P=P,
                 logFC_value=logFC_value,P_value=P_value,
                 feature_type=feature_type,colors=colors)
  save(DEG_list,file = "data/DEG_list.rda")
  write.csv(DEG,'Table/DEG.csv')
  
  exp_plot=log2(cpm(exp)+1)#用于后续绘图
  }
}else{
  ####limma####
  library(limma)
  library(dplyr)
  group_deg <- factor(group,levels = c(control,case))
  ## 实验设计矩阵
  design <- model.matrix(~ 0 + group_deg)
  rownames(design) <- colnames(exp)
  colnames(design) <- levels(group_deg)
  ## 线性建模
  fit <- lmFit(exp,design)
  cont.matrix <- makeContrasts(contrasts = paste0(case,'-',control), levels = design)#设置实验组和对照组contrasts = c("实验组-对照组")
  fit2 <- contrasts.fit(fit, cont.matrix)
  ## 经验贝叶斯调整
  fit2 <- eBayes(fit2)
  ## 筛选差异基因
  dif <- topTable(fit2, coef = 1, n = Inf)
  dif <- na.omit(dif)
  dif$gene_symbol <- row.names(dif)
  if(P=="p.adj"){
    dif.up <- dif %>%
      dplyr::filter(logFC > logFC_value & adj.P.Val < P_value)#筛选上调基因，logFC可调整，P值不能大于0.05
    dif.down <- dif %>%
      dplyr::filter(logFC < (-logFC_value) & adj.P.Val < P_value)#筛选下调基因，logFC可调整，P值不能大于0.05
  }else{
    dif.up <- dif %>%
      dplyr::filter(logFC > logFC_value & P.Value < P_value)#筛选上调基因，logFC可调整，P值不能大于0.05
    dif.down <- dif %>%
      dplyr::filter(logFC < (-logFC_value) & P.Value < P_value)#筛选下调基因，logFC可调整，P值不能大于0.05
  }
  #添加差异基因信息
  dif <- dif %>% dplyr::mutate(group = dplyr::case_when(
    gene_symbol %in% dif.up$gene_symbol ~ "up",
    gene_symbol %in% dif.down$gene_symbol ~ "down",
    TRUE ~ "no"
  ))
  
  DEG=data.frame(name=dif$gene_symbol,
                 logFC=dif$logFC,
                 P.value=dif$P.Value,
                 P.adj=dif$adj.P.Val,
                 group=dif$group)
  ####输出数据####
  dir.create("data")
  dir.create("Table")
  print(table(DEG$group))
  
  DEG_list<-list(method="limma",
                 case=case,
                 control=control,
                 DEG=DEG,P=P,
                 logFC_value=logFC_value,P_value=P_value,
                 feature_type=feature_type,colors=colors)
  save(DEG_list,file = "data/DEG_list.rda")
  write.csv(DEG,'Table/DEG.csv')
  
  exp_plot=exp#用于后续绘图
}
####绘图####
dir.create("plot")
top10_padj = DEG[DEG$group!='no',] %>% group_by(group) %>% slice_max(n = 5,order_by = -log10(P.adj),with_ties = F)#选择P值最小的差异基因上下调各五个
top10_p = DEG[DEG$group!='no',] %>% group_by(group) %>% slice_max(n=5,order_by = -log10(P.value),with_ties = F)
####火山图####
library(ggplot2)
if(P=="p")label=top10_p else label=top10_padj#选择P值最小的差异基因上下调各五个
DEG$Label=ifelse(DEG$name%in%label$name,DEG$name,NA)#添加一列存放选取出来的差异基因
Significant=factor(DEG$group,levels = c('up','down',"no"))
dp<-ifelse(P=='p',3,4)
library(ggplot2)
library(ggrepel)
p1<-ggplot(DEG, aes(logFC, -log10(DEG[,dp])))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c(colors[2],colors[1],"grey"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-logFC_value,logFC_value), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(plot.title = element_blank(),plot.subtitle = element_blank(),
        plot.background = element_blank(),plot.margin = margin(t=1,r=1,b=1,l=1,unit="pt"),
        panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6,colour = "black"),axis.title=element_text(size=7),
        legend.background = element_blank(),legend.key = element_blank(),
        legend.title = element_text(size = 7),legend.text = element_text(size = 7),
        legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),legend.box.spacing = unit(3,"pt"))+
  str(DEG, max.level = c(-1, 1))+geom_text_repel(aes(label=Label),size=2.5)#添加注释文本

while (!is.null(dev.list()))  dev.off()
pdf(paste0('plot/',case,'_vs_',control,'_volplot.pdf'),width = 8/2.54,height = 8/2.54)
print(p1)
dev.off()
####热图1####
if(P=="p")top10=top10_p else top10=top10_padj
group_annotate=group_deg[order(group_deg)]
annotate<-as.data.frame(group_annotate)
mat<-exp_plot[top10$name,rownames(annotate)]
colnames(annotate)="sample"

annotate_color<-list(sample=c(colors[1],colors[2]))
names(annotate_color$sample)<-levels(group_annotate)

mat=t(scale(t(mat)))
mat[mat > 3] = 3 #限定上限，使表达量大于3的等于3
mat[mat < -3] = -3 #限定下限，使表达量小于-3的等于-3

p<-pheatmap1(mat,annotation_col  = annotate,scale = 'none',
            show_rownames = T,annotation_colors = annotate_color,
            show_colnames = F,cluster_cols = F,cluster_rows = T,fontsize = 5,name = " ",
            color = colorRampPalette(colors = c(colors[1],"white",colors[2]))(50))

while (!is.null(dev.list()))  dev.off()
pdf(paste0('plot/',case,'-',control,'_top10_heatmap.pdf'),width = (max(nchar(rownames(mat)))/8+8)/2.54,height = 8/2.54)
print(p)
dev.off()

####热图2####

mat<-exp_plot[DEG$name[DEG$group!='no'],rownames(annotate)]
mat=t(scale(t(mat)))
mat[mat > 3] = 3 #限定上限，使表达量大于3的等于3
mat[mat < -3] = -3 #限定下限，使表达量小于-3的等于-3

p<-pheatmap1(mat,annotation_col  = annotate,scale = 'none',
             show_rownames = F,annotation_colors = annotate_color,border_color = NA,
             show_colnames = F,cluster_cols = F,cluster_rows = T,fontsize = 5,name = " ",
             color = colorRampPalette(colors = c(colors[1],"white",colors[2]))(50))

while (!is.null(dev.list()))  dev.off()
pdf(paste0('plot/',case,'-',control,'_allDEG_heatmap.pdf'),width = 8/2.54,height = 8/2.54)
print(p)
dev.off()

####箱线图####
p<-plot_boxplot(data = as.data.frame(t(exp_plot)),
                group = group,
                feature = top10$name,
                colors = colors,#颜色参数
                flip = flip,#是否进行翻转
                lab = lab,#以x，y的顺序输入坐标轴标题
                feature_type=feature_type)#feature的类型
pdf(paste0('plot/',case,'-',control,'_top10_boxplot.pdf'),width = 16/2.54,height = (max(nchar(rownames(exp_plot)))/8+8)/2.54)
print(p)
dev.off()
####小提琴图####
p<-plot_vioplot(data = as.data.frame(t(exp_plot)),
                group = group,
                feature = top10$name,
                colors = colors,#颜色参数
                flip = flip,#是否进行翻转
                box = F,
                lab = lab,#以x，y的顺序输入坐标轴标题
                feature_type=feature_type)#feature的类型
pdf(paste0('plot/',case,'-',control,'_top10_vioplot.pdf'),width = 16/2.54,height = (max(nchar(rownames(exp_plot)))/8+8)/2.54)
print(p)
dev.off()

setwd("../")
