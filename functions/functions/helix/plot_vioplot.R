####小提琴图####
library(ggplot2)
library(dplyr)
library(ggpubr)
plot_vioplot<-function(data=NULL,feature=NULL,group=NULL,colors=NULL,flip=F,box=F,feature_type=NULL,lab=NULL){
  #数据检测####
  if(is.null(data))stop("缺少data")
  if(is.null(group))stop("缺少group")
  if(is.null(colors))stop("缺少colors")
  if(is.null(feature_type))stop("缺少feature_type")
  if(all(!is.null(lab),length(lab)<2,length(lab)>2))stop("请以x、y的顺序设置标题名称")
  if(!all(rownames(data)==names(group)))stop("data和group样本名未完全对应")
  #数据整理####
  if(!is.null(feature)){
    if(length(feature)==1){
      data=data.frame(data[,feature],row.names = rownames(data))
      colnames(data)=feature
    }else data=data[,feature]
  }
  plot_data=data%>%rownames_to_column("Sample")%>%
    gather(key = key,value = value,-Sample)
  plot_data$group=group[plot_data$Sample]
  #秩和检验####
  group_class<-unique(group)
  data_1=as.data.frame(data[names(group)[group==group_class[1]],])
  colnames(data_1)=colnames(data)
  data_2=as.data.frame(data[names(group)[group==group_class[2]],])
  colnames(data_2)=colnames(data)
  for(i in colnames(data)){
    res1<-wilcox.test(data_1[,i],data_2[,i])
    if(i == colnames(data)[1]){
      wilcoxtest<-data.frame(feature=i,p.value=res1$p.value,up_group=ifelse(res1$p.value>0.05,"no sig",ifelse(median(data_1[,i])>median(data_2[,i]),group_class[1],group_class[2]))) 
    }else{
      tmp<-data.frame(feature=i,p.value=res1$p.value,up_group=ifelse(res1$p.value>0.05,"no sig",ifelse(median(data_1[,i])>median(data_2[,i]),group_class[1],group_class[2])))
      wilcoxtest<-rbind(wilcoxtest,tmp)
    }
  }
  vioplot_text=paste0("如图（Figure）所示，共有",sum(wilcoxtest$p.value<0.05&wilcoxtest$up_group==group_class[1]),
                      "个",feature_type,"在",group_class[1],"组中显著升高，",sum(wilcoxtest$p.value<0.05&wilcoxtest$up_group==group_class[2]),
                      "个",feature_type,"在",group_class[2],"组中显著升高（Table）")
  #输出文字和表格####
  vioplot_test_list<<-list(vioplot_text=vioplot_text,
                           vioplot_wilcoxtest=wilcoxtest)
  #绘图####
  p<-ggplot(plot_data,aes(key,value,fill = group))+
    geom_violin(alpha = 0.4,trim=F)+ 
    scale_fill_manual(values = c(colors[1],colors[2]))+
    stat_compare_means(aes(group = group,label = ..p.signif..),bracket.size = 0.6, size = 3,
                       label.y = max(plot_data$value)*0.9+0.1*min(plot_data$value), #p值位置
                       hide.ns = T,#隐藏ns
                       method = ifelse(length(unique(plot_data$group))==2,"wilcox.test","kruskal.test"))+
    theme(plot.title = element_blank(),plot.subtitle = element_blank(),
          plot.background = element_blank(),plot.margin = margin(t=2,r=2,b=2,l=2,unit="pt"),
          panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1,size = 7,colour = "black"),
          axis.title = element_text(size = 7,colour = "black"),
          axis.text.y = element_text(size = 6,colour = "black"),
          legend.background = element_blank(),legend.key = element_blank(),
          legend.title = element_text(size = 7),legend.text = element_text(size = 7),
          legend.position = "top",legend.box.spacing = unit(3,"pt"))
  
  #其他修饰####
  #旋转箱线图
  if(flip){
    p=p+coord_flip() 
  }
  #增加箱线图
  if(box){
    p=p+geom_boxplot(width=0.1,position = position_dodge(0.9),alpha = 0.4,outlier.size = 0)
  }
  if(!is.null(lab)){
    p=p+labs(x = lab[1], y = lab[2])
  }else p=p+labs(x = "", y = "Expression Level")
  #输出图####
  return(p)
}

