####泛癌分析####
####数据检测和处理####
dir.create('Pan_cancer')
dir.create("Pan_cancer/plot")
dir.create("Pan_cancer/Table")
dir.create("Pan_cancer/data")
####关键基因泛癌箱线图####

for(gene in feature){
  
  dir.create(paste0('Pan_cancer/plot/',gene))
  
  box_data<-data.frame()
  for(i in names(TCGA_list)){
    data<-data.frame(
      sample=TCGA_list[[i]]$group$sample,
      value=t(TCGA_list[[i]]$exp[gene,TCGA_list[[i]]$group$sample])[,1],
      group=TCGA_list[[i]]$group$group,
      type=i
    )
    box_data<-rbind(box_data,data)
  }
  
  library(ggplot2)
  library(ggpubr)
  boxplot_pancancer<-ggplot(data = box_data,aes(x=type,y=value,fill=group))+
    geom_boxplot(outlier.shape = 21,outlier.color = NA,outlier.size = 0.8,
                 #notch=TRUE,notchwidth=0.8, #在箱子上生成槽口,notchwidth越小则越往里凹
                 alpha = 0.4)+  
    labs(y=gene,x='Type')+
    scale_fill_manual(values = c(colors[2],colors[1]))+
    stat_compare_means(aes(group = group,label = ..p.signif..),bracket.size = 0.6, size = 2.5,
                       label.y = max(box_data$value)*0.9+min(box_data$value)*0.1, #p值位置
                       hide.ns = T,#隐藏ns
                       method = "wilcox.test")+
    theme(plot.title = element_blank(),plot.subtitle = element_blank(),
          plot.background = element_blank(),
          panel.border = element_rect(size = 0.5,fill = NA,colour = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1,size = 7,colour = "black"),
          axis.title = element_text(size = 7,colour = "black"),
          axis.text.y = element_text(size = 6,colour = "black"),
          legend.background = element_blank(),legend.key = element_blank(),
          legend.title = element_text(size = 7),legend.text = element_text(size = 7),
          legend.position = "top",
          legend.box.spacing = unit(3,"pt"))#+coord_flip() #旋转箱线图
  
  pdf(paste0('Pan_cancer/plot/',gene,'/',gene,'_pancancer_boxplot.pdf'),width = 17/2.54,height = 10/2.54)
  print(boxplot_pancancer)
  dev.off()
  
  ####单因素####
  
  library(survival)
  cox_data<-data.frame()
  for(i in names(TCGA_list)){
    data<-cbind(t(TCGA_list[[i]]$exp_tumor[gene,]),TCGA_list[[i]]$clinical[,c('fustat','futime')])
    Bcox<-coxph(Surv(futime, fustat)~as.numeric(as.character(data[,gene]))>median(as.numeric(as.character(data[,gene]))),data=data)
    summcph<-summary(Bcox)
    sur.res<-data.frame(
      "Cancer"=i,
      "HR"=round(summcph$conf.int[1],2),
      "CI5"=round(summcph$conf.int[3],2),
      "CI95"=round(summcph$conf.int[4],2),
      "P"=ifelse(summcph$coefficients[,5] < 0.05,"<0.05",round(summcph$coefficients[,5],3))
    )
    
    rownames(sur.res)<-i
    cox_data<-rbind(cox_data,sur.res)
  }
  
  #将95%CI连接起来
  cox_data$HR.CI95 <- paste0(cox_data$HR," (",cox_data$CI5,'-',cox_data$CI95,")")
  result <- cox_data[,c(1:4,6,5)]
  
  #给参考变量插入空行
  ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
  #插入空行，形成一个新表
  result1<-data.frame("Cancer", NA, NA, NA, "HR(95%CI)","P.value")
  colnames(result1)<-colnames(result)
  
  result1<-rbind(result1,
                 result,
                 c(NA, NA, NA, NA, NA,NA))
  
  rownames(result1)<-1:nrow(result1)
  
  ####图形进行美化
  library(forestplot)
  hrzl_lines<-list(gpar(lty=1,lwd=2),#表头上方添加实线
                   gpar(lty=2),#表头下方添加虚线
                   gpar(lwd=2,lty=1,columns=c(1:4)))#最后一行下方添加实线
  
  names(hrzl_lines)<-c('1','2',nrow(result1)+1)
  
  is_summary<-is.na(result1[,2]) #按顺序指定每行是否加粗，T为加粗，F不加粗
  is_summary[1]=TRUE
  
  fig2<- forestplot(result1[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                    mean=as.numeric(result1[,2]),   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                    lower=as.numeric(result1[,3]),  #告诉函数表格第3列为5%CI，
                    upper=as.numeric(result1[,4]),  #表格第4列为95%CI，它俩要化作线段，穿过方块
                    zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                    boxsize=0.5,       #设置小黑块的大小
                    graph.pos= "right" ,
                    hrzl_lines=hrzl_lines, #最后一行下方添加实线
                    graphwidth = unit(.25,"npc"),
                    #xticks=c(-1,1,2,3,4,5) , #森林图刻度
                    is.summary=is_summary, #按顺序指定每行是否加粗，T为加粗，F不加粗
                    txt_gp=fpTxtGp(
                      label=gpar(cex=0.5),
                      ticks=gpar(cex=0.5), 
                      xlab=gpar(cex=0.5), 
                      title=gpar(cex=0.5)),
                    lwd.zero=1,
                    lwd.ci=1.5,
                    lwd.xaxis=2, 
                    lty.ci=1.5,
                    ci.vertices =T,
                    ci.vertices.height=0.2, 
                    clip=c(0,5),
                    ineheight=unit(9, 'mm'), 
                    line.margin=unit(9, 'mm'),
                    colgap=unit(2, 'mm'),
                    fn.ci_norm="fpDrawDiamondCI", 
                    title=gene,
                    col=fpColors(box =colors, 
                                 lines =colors, 
                                 zero = "black"))
  while (!is.null(dev.list()))  dev.off()
  pdf(paste0("Pan_cancer/plot/",gene,'/',gene,"_单因素Cox回归森林图.pdf"),height = 21/2.54,width = 12/2.54)
  print(fig2)
  dev.off()
  
  ####肿瘤正常对应图####
  # tumor<-"TCGA-STAD"
  # box_data1<-data.frame(sample=TCGA_list[[tumor]]$sample_pa$sample,
  #                       sample1=TCGA_list[[tumor]]$sample_pa$sample1,
  #                       group=TCGA_list[[tumor]]$sample_pa$group,
  #                       value=t(TCGA_list[[tumor]]$exp[gene,TCGA_list[[tumor]]$sample_pa$sample])[,1])
  # 
  # ggplot(data = box_data1,aes(x=group,y=value,fill=group))+
  #   geom_point(aes(color=group),size=3)+
  #   geom_line(aes(group=sample1),lwd=1)+
  #   labs(y=gene,x='Type')+
  #   scale_fill_manual(values = c(colors[2],colors[1]))+
  #   scale_color_manual(values = c(colors[2],colors[1]))+
  #   stat_compare_means(aes(group = group,label = ..p.signif..),bracket.size = 0.6, size = 3,
  #                      label.y = max(box_data$value)-0.02, #p值位置
  #                      hide.ns = T,#隐藏ns
  #                      method = "wilcox.test")+
  #   theme(plot.title = element_blank(),plot.subtitle = element_blank(),
  #         plot.background = element_blank(),plot.margin = margin(t=2,r=2,b=2,l=2,unit="pt"),
  #         panel.border = element_rect(size = 0.5,fill = NA,colour = "black"),
  #         panel.background = element_blank(),
  #         panel.grid = element_blank(),
  #         axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1,size = 7,colour = "black"),
  #         axis.title = element_text(size = 7,colour = "black"),
  #         axis.text.y = element_text(size = 6,colour = "black"),
  #         legend.background = element_blank(),legend.key = element_blank(),
  #         legend.title = element_text(size = 7),legend.text = element_text(size = 7),
  #         legend.position = "top",legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),
  #         legend.box.spacing = unit(3,"pt"))#+coord_flip() #旋转箱线图
  # 
  
  #带显著性的相关性热图
  library(tidyverse)
  library(reshape2)
  library(Hmisc)
  
  for(i in 1:length(TCGA_list)){
    exp<-t(TCGA_list[[i]]$ssGSEA)
    table(rownames(exp)==colnames(TCGA_list[[i]]$exp_tumor))
    exp=cbind(exp,t(TCGA_list[[i]]$exp_tumor)[,gene])
    colnames(exp)[ncol(exp)]<-names(TCGA_list)[i]
    
    rcorr_list<-rcorr(exp,type = "pearson")#"pearson","spearman"
    cor_ma<-rcorr_list$r
    cor_ma1<-data.frame(cor_ma[-29,29])
    colnames(cor_ma1)<-str_split(names(TCGA_list)[i],pattern = '-',simplify = T)[,2]
    if(i==1)cor_ma2<-cor_ma1 else cor_ma2<-cbind(cor_ma2,cor_ma1)
    
    p_ma<-rcorr_list$P
    p_ma1<-data.frame(p_ma[-29,29])
    colnames(p_ma1)<-str_split(names(TCGA_list)[i],pattern = '-',simplify = T)[,2]
    if(i==1)p_ma2<-p_ma1 else p_ma2<-cbind(p_ma2,p_ma1)
  }
  
  data_cor <- t(cor_ma2) %>% as.data.frame() %>%
    rownames_to_column("x") %>% 
    gather(key = y,value = cor,-x)
  
  data_p <- t(p_ma2) %>% as.data.frame() %>%
    rownames_to_column("x") %>% 
    gather(key = y,value = p,-x)
  
  data_all<-merge(data_cor,data_p,by=c('x','y'))
  
  data_all$p_text<-ifelse(data_all$p<0.0001,"****",
                          ifelse(data_all$p<0.001,"***",
                                 ifelse(data_all$p<0.01,"**",
                                        ifelse(data_all$p<0.05,"*",""))))
  
  p<-ggplot(data_all, aes(x, y)) + 
    geom_tile(aes(fill = cor), colour = "white", size = 0.1)+
    geom_text(aes(label=p_text),col ="black",size = 2) +
    scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") + 
    scale_color_gradient2(low = "#5C5DAF",mid = "#FCB886",high = "#EA2E2D") + 
    theme_minimal() + # 不要背景
    theme(#axis.title.x=element_blank(), # 去掉 title
          axis.ticks.x=element_blank(), # 去掉x 轴
          axis.title.y=element_blank(), # 去掉 y 轴
          axis.text.x = element_text(size = 6, face = "plain",angle = 90), # 调整x轴文字
          axis.text.y = element_text(size = 6, face = "plain"),#调整y轴文字
          legend.text = element_text(size = 7, face = "plain"),
          legend.title = element_text(size = 7, face = "plain"),
          panel.grid.major =  element_blank(),
          legend.position = 'right') + 
    labs(fill =paste0("* p < 0.05","\n\n","** p <0.01","\n\n","*** p <0.001","\n\n","**** p <0.0001","\n\n","Correlation"),
         x=paste0(gene," Expression"))# 修改 legend 内容
  while (!is.null(dev.list()))  dev.off()
  pdf(paste0('Pan_cancer/plot/',gene,'/',gene,'_immune_cells_cor_heatmap.pdf'),width =8.5 ,height = 6)
  print(p)
  dev.off()
  
}
