####基因表达情况验证####
dir.create("Gene_expression")
dir.create("Gene_expression/data")
dir.create("Gene_expression/plot")
dir.create("Gene_expression/Table")
#数据检测和处理####
data=data[,feature]
group_annotate=factor(group,levels = c(control,case))

#热图####
group_annotate=group_annotate[order(group_annotate)]
annotate<-as.data.frame(group_annotate)
mat<-data[rownames(annotate),]
colnames(annotate)="sample"

annotate_color<-list(sample=c(colors[1],colors[2]))
names(annotate_color$sample)<-levels(group_annotate)

mat=t(scale(mat))
mat[mat > 3] = 3 #限定上限，使表达量大于3的等于3
mat[mat < -3] = -3 #限定下限，使表达量小于-3的等于-3

p<-pheatmap1(mat,annotation_col  = annotate,scale = 'none',
             show_rownames = T,annotation_colors = annotate_color,
             show_colnames = F,cluster_cols = F,cluster_rows = T,fontsize = 5,name = " ",
             color = colorRampPalette(colors = c(colors[1],"white",colors[2]))(50))

while (!is.null(dev.list()))  dev.off()
pdf(paste0('Gene_expression/plot/',case,'-',control,'_heatmap.pdf'),width = (max(nchar(rownames(mat)))/8+8)/2.54,height = 8/2.54)
print(p)
dev.off()
#箱线图####
p<-plot_boxplot(data = data,
                group = group,
                feature = feature,
                colors = colors,#颜色参数
                flip = flip,#是否进行翻转
                lab = lab,#以x，y的顺序输入坐标轴标题
                feature_type=feature_type)#feature的类型

while (!is.null(dev.list()))  dev.off()
pdf(paste0('Gene_expression/plot/',case,'-',control,'_boxplot.pdf'),width = (max(nchar(rownames(mat)))/8+8)/2.54,height = 8/2.54)
print(p)
dev.off()

#小提琴图####
p<-plot_vioplot(data = data,
                group = group,
                feature = feature,
                colors = colors,#颜色参数
                flip = flip,#是否进行翻转
                box = T,
                lab = lab,#以x，y的顺序输入坐标轴标题
                feature_type=feature_type)#feature的类型

while (!is.null(dev.list()))  dev.off()
pdf(paste0('Gene_expression/plot/',case,'-',control,'_vioplot.pdf'),width = (max(nchar(rownames(mat)))/8+8)/2.54,height = 8/2.54)
print(p)
dev.off()

#ROC曲线####
data1<-data
data1$group=group
library(pROC)
library(ggplot2)

dir.create("Gene_expression/plot/ROC")

#data为包涵分类信息（type）和需要研究的基因的表达矩阵
#type为要研究的分类（数值型的二分类，如1为肿瘤0为正常）
#gene为要研究的基因
auc_table=data.frame()
ROC_list<-list()
hubgene<-feature
for(i in hubgene){
  rocobj1 <- roc(data1$group, data1[,i],    #list of 15
                 smooth = F)       # 曲线是否光滑，当光滑时，无法计算置信区间 
  # 计算AUC值
  auc<-auc(rocobj1)[1]
  auc_table=rbind(auc_table,data.frame(Gene=i,AUC=auc))
  if (auc>0.6){
    auc_text<-paste0("AUC = ",round(auc,4))
    # 绘图
    p_roc=ggroc(rocobj1,
                color=colors[3],
                size=1,
                legacy.axes = F # FALSE时 横坐标为1-0 specificity；TRUE时 横坐标为0-1 1-specificity
    )+
      theme_classic()+
      geom_ribbon(aes(x=specificity, ymin = 0, ymax = sensitivity),
                  fill = colors[3], alpha = .6)+
      # geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),        # 绘制对角线
      #              colour='grey', 
      #              linetype = 'dotdash'
      # ) +
      annotate('text',x = 0.25,y = 0.25,label=paste0(i,'\n',auc_text),size=2.5)+
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
    pdf(paste0('Gene_expression/plot/ROC/',i,"_ROC.pdf"),width = 6/2.54,height = 6/2.54)
    print(p_roc)
    while (!is.null(dev.list()))  dev.off()
    ROC_list[[i]]=p_roc}
}  
write.csv(auc_table,"Gene_expression/Table/hubgene_AUC.csv")
#devtools::install_github("thomasp85/patchwork@v1.1.0")
library(patchwork)
library(ggplotify)

hubgene=auc_table$Gene[auc_table$AUC>0.6]
j_roc=0
for (i in 1:length(hubgene)) {
  if (i%%12==1) {
    patchwork=paste0('ROC_list[[hubgene[',i,']]]+')
    j_roc=j_roc+1
  }else if (i==length(hubgene)|i%%12==0) {
    if (i%%12<3&i%%12!=0) {
      patchwork=paste0(patchwork,'ROC_list[[hubgene[',i,']]]')
      eval(parse(text = paste0('patchwork=',patchwork)))
      p_roc_p<-patchwork + plot_layout(ncol = i%%12)+
        plot_annotation(tag_levels = 'A') & 
        theme(plot.tag = element_text(size = 15))
    }else{
      patchwork=paste0(patchwork,'ROC_list[[hubgene[',i,']]]')
      eval(parse(text = paste0('patchwork=',patchwork)))
      p_roc_p<-patchwork + plot_layout(ncol = 3)+
        plot_annotation(tag_levels = 'A') & 
        theme(plot.tag = element_text(size = 15))
    }
    
    height_roc=21
    width_roc=17
    if(i==length(hubgene)){height_roc=ifelse(length(hubgene)%%12%%3==0,length(hubgene)%%12/3/4*height_roc,
                                             (floor(length(hubgene)%%12/3)+1)/4*height_roc)
    if(length(hubgene)%%12==0)height_roc=21
    width_roc=ifelse(length(hubgene)%%12<3&length(hubgene)%%12!=0,length(hubgene)%%12/3*width_roc,width_roc)
    }
    pdf(paste0('Gene_expression/plot/ROC_plot_',j_roc,'.pdf'),width = width_roc/2.54,height = height_roc/2.54)
    print(p_roc_p)
    dev.off()
  }else{patchwork=paste0(patchwork,'ROC_list[[hubgene[',i,']]]+')}
  
}
