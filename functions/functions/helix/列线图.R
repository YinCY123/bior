
####nomogram####
#参数检测
if(!is.numeric(data[,event_use])){
  if(is.factor(data[,event_use])){
    if(length(levels(data[,event_use]))!=2)stop("结局数据分组数不对，请检查")
    if(min(as.numeric(data[,event_use]))>0)data[,event_use]=as.numeric(data[,event_use])-1 else 
      data[,event_use]=as.numeric(data[,event_use])
  } else stop("结局列event_use应为0，1的数值格式或因子格式，且对照组应为数值小的项")
}

#函数定义
reduce_vector_uniformly <- function(original_vector, new_length) {  
  if (length(original_vector) <= new_length) {  
    return(original_vector)  
  }  
  
  # 保留首尾值  
  new_vector <- c(original_vector[1], rep(NA, new_length - 2), original_vector[length(original_vector)])  
  
  # 计算理想的等分点位置  
  intervals <- seq(from = 0, to = 1, length.out = new_length)  
  intervals <- intervals[-c(1, length(intervals))]  # 移除首尾点，因为已经在新向量中了  
  
  # 在原向量中找到与等分点最接近的值  
  for (i in 1:length(intervals)) {  
    target <- intervals[i] * (length(original_vector) - 1) + 1  # 将[0,1]区间映射到原向量的索引  
    indices <- 1:length(original_vector)  
    closest_index <- which.min(abs(indices - target))  # 找到最接近的索引  
    new_vector[i + 1] <- original_vector[closest_index]  # 将找到的值放入新向量  
  }  
  
  return(new_vector)  
}  

dir.create('Nomogram')
dir.create("Nomogram/data")
dir.create("Nomogram/plot")
dir.create("Nomogram/Table")
Nomogram_list=list()

if(!diagnosis){
  ####预后模型####
  y<- Surv(time = data[,time_use],event = data[,event_use])#1为感兴趣事件
  FML=as.formula(paste0("y~",paste(feature,collapse = "+")))
  
  dd=datadist(data)#data必须是data.frame
  options(datadist="dd")
  fit<-cph(FML,data=data,x=TRUE,y=TRUE,surv=TRUE)
  mode(data)
  
  survival<-Survival(fit)
  survival1<-function(x)survival(years[1],x)
  survival2<-function(x)survival(years[2],x)
  survival3<-function(x)survival(years[3],x)
  
  nomo<-nomogram(fit,fun=list(survival1,survival2,survival3),
                 funlabel =paste0(years," Year Survival"))
  
  #C-index检验预测效果
  B<-rcorrcens(y~ predict(fit),data=data)
  C_Index=(1-B[1])
  #调整刻度数目
  
  for(i in names(nomo)[grepl("Survival",names(nomo))]){
    for(j in 1:3){
      n=length(nomo[[i]][[j]])
      nomo[[i]][[j]]=nomo[[i]][[j]][reduce_vector_uniformly(1:n,nomo_n)]
    }
  }
  
  Nomogram_list[["nomogram"]]=nomo
  Nomogram_list[["C_Index"]]=C_Index
  
  while (!is.null(dev.list()))  dev.off()
  pdf("Nomogram/plot/nomplot.pdf",width=6,height=6)
  par(family="sans",ps="6") #sans为Arial字体，serif为新罗马字体 ,ps为字体磅值
  plot(nomo, lplabel="Linear Predictor",
       xfrac=.25, #左边标题与右边图形间隔
       label.every = 1,col.grid = gray(c(0.8, 0.95)), #对应上方points的线条
       cex.var = 1.2,cex.axis = 1,cex.lab = 1.2, #cex文本属性
       lty=1, #lty指定线型
       lwd=5 #lwd改变线条粗细
  )
  title(main="")
  dev.off()
  
  ####画校准曲线calibration curve
  fit1<-cph(FML,data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[1])
  fit2<-cph(FML,data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[2])
  fit3<-cph(FML,data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[3])
  
  cal1<-calibrate(fit1, cmethod="KM",method="boot", u=years[1], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关
  
  cal2<-calibrate(fit2, cmethod="KM",method="boot", u=years[2], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关
  
  cal3<-calibrate(fit3, cmethod="KM",method="boot", u=years[3], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关
  
  while (!is.null(dev.list()))  dev.off()
  pdf("Nomogram/plot/Calibration.pdf",width=4,height=4.5)
  par(family="sans",ps="6",mar=c(6,4,1,1),mgp=c(1.5, 0.5, 0),tck = -0.02)
  plot(cal1,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[1],col=colors[1],lwd=1,pch=NA,par.corrected = list(col=colors[1],lwd=1,pch=15),
       xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")
  
  par(new=TRUE,family="sans",ps="0")
  plot(cal2,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[2],col=colors[2],lwd=1,axes=FALSE,pch=NA,par.corrected = list(col=colors[2],lwd=1,pch=15),
       xlab="",ylab="", add=T)
  
  par(new=TRUE,family="sans",ps="6")
  plot(cal3,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[3],col=colors[3],lwd=1,axes=FALSE,pch=NA,par.corrected = list(col=colors[3],lwd=1,pch=15),
       xlab="",ylab="", add=T)  #errbar.col定义误差线的颜色，col定义校准曲线的颜色
  
  legend("topleft", inset=0.05,paste0(years,'-year'), col=colors[1:3], lty =1,lwd =1, bty = "n")
  
  dev.off()
  
  ####用ROC曲线验证
  library(rmda)
  
  nomoRisk=predict(fit, data=data, type="lp")
  data$nomoRisk=nomoRisk
  
  rt <- data
  ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
                 marker=rt$nomoRisk, cause=1,
                 weighting='aalen',
                 times=years, ROC=TRUE)
  while (!is.null(dev.list()))  dev.off()
  pdf(file="Nomogram/plot/ROC.pdf",width=5,height=5)
  plot(ROC_rt,time=years[1],col=colors[3],title=FALSE,lwd=2)
  plot(ROC_rt,time=years[2],col=colors[2],add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=years[3],col=colors[1],add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at ',years[1],' years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at ',years[2],' years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at ',years[3],' years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c(colors[3],colors[2],colors[1]),lwd=2,bty = 'n')
  dev.off()
  Nomogram_list[["data"]]=data
  save(Nomogram_list,file = "Nomogram/data/Nomogram_list.rda")
}else{
  ####诊断模型####
  
  library(rms)
  y<- data[,event_use]#1为感兴趣事件
  FML=as.formula(paste0("y~",paste(feature,collapse = "+")))
  
  fit<-lrm(FML,data = data,x=T,y=T)
  
  dd<-datadist(data)
  options(datadist='dd')
  nomo<-nomogram(fit,fun = plogis,conf.int = F,
                 funlabel = "Risk of Disease", lp=F)
  
  for(i in names(nomo)[grepl("Disease",names(nomo))]){
    for(j in 1:3){
      n=length(nomo[[i]][[j]])
      nomo[[i]][[j]]=nomo[[i]][[j]][reduce_vector_uniformly(1:n,nomo_n)]
    }
  }
  
  pdf("Nomogram/plot/nomplot.pdf",width=8,height=6)
  par(family="sans",ps="6") #sans为Arial字体，serif为新罗马字体 ,ps为字体磅值
  plot(nomo, total.points.label='Total Points',
       xfrac=.25, #左边标题与右边图形间隔
       label.every = 1,col.grid = gray(c(0.8, 0.95)), #对应上方points的线条
       cex.var = 1.2,cex.axis = 1,cex.lab = 1.2, #cex文本属性
       lty=1, #lty指定线型
       lwd=5 #lwd改变线条粗细
  )
  dev.off()
  
  Nomogram_list[["nomogram"]]=nomo
  Nomogram_list[["fit"]]=fit
  P1 <- predict(fit,type = 'fitted')  ##获得预测概率值
  
  calibrate_table=val.prob(P1,y) 
  write.csv(calibrate_table,"Nomogram/Table/calibrate_table.csv")
  
  cal<-calibrate(fit,method = 'boot',B=500)
  while (!is.null(dev.list()))  dev.off()
  pdf('Nomogram/plot/Calibration.pdf',width = 7,height = 7)
  plot(cal,
       xlim = c(0,1),
       xlab = "Predicted Probability",
       ylab = "Observed Probability",
       legend = FALSE,
       subtitles = FALSE)
  abline(0,1,col = "black",lty = 2,lwd = 2)
  lines(cal[,c("predy","calibrated.orig")], type = "l",lwd = 2,col=colors[1],pch =16)
  lines(cal[,c("predy","calibrated.corrected")], type = "l",lwd = 2,col=colors[2],pch =16)
  legend(0.55,0.35,
         c("Apparent","Ideal","Bias-corrected"),
         lty = c(2,1,1),
         lwd = c(2,1,1),
         col = c("black",colors[1],colors[2]),
         bty = "n") # "o"为加边框
  dev.off()
  
  Nomogram_list[["calibrate"]]=cal
  Nomogram_list[["calibrate_table"]]=calibrate_table
  #单基因ROC
  data1=data
  data1$nomoRisk=P1
  
  library(pROC)
  library(ggplot2)
  
  #列线图ROC
  rocobj1 <- roc(data1[,event_use], data1[,"nomoRisk"],    #list of 15
                 smooth = F)       # 曲线是否光滑，当光滑时，无法计算置信区间 
  # 计算AUC值
  auc<-auc(rocobj1)[1]
  auc_table=data.frame(Gene="nomoRisk",AUC=auc)
  auc_text<-paste0("AUC = ",round(auc,4))
  # 绘图
  p_roc=ggroc(rocobj1,
              color=colors[1],
              size=1,
              legacy.axes = F # FALSE时 横坐标为1-0 specificity；TRUE时 横坐标为0-1 1-specificity
  )+
    theme_classic()+
    geom_ribbon(aes(x=specificity, ymin = 0, ymax = sensitivity),
                fill = colors[1], alpha = .6)+
    # geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),        # 绘制对角线
    #              colour='grey', 
    #              linetype = 'dotdash'
    # ) +
    annotate('text',x = 0.25,y = 0.25,label=paste0("nomoRisk",'\n',auc_text),size=2.5)+
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  pdf(paste0("Nomogram/plot/ROC_use.pdf"),width = 5.1653/2.54,height = 4.5991/2.54)
  print(p_roc)
  while (!is.null(dev.list()))  dev.off()
  pdf(paste0("Nomogram/plot/ROC.pdf"),width = 6/2.54,height = 6/2.54)
  print(p_roc)
  while (!is.null(dev.list()))  dev.off()
  Nomogram_list[["data"]]=data1
  save(Nomogram_list,file = "Nomogram/data/Nomogram_list.rda")
}

Nomogram_valid<-function(data=NULL,feature=feature){
  load("Nomogram/data/Nomogram_list.rda")
  fit=Nomogram_list$fit
  P1=predict(fit, newdata = data,type = 'fitted')
  #单基因ROC
  data1=data
  data1$nomoRisk=P1
  
  library(pROC)
  library(ggplot2)
  
  #列线图ROC
  rocobj1 <- roc(data1[,event_use], data1[,"nomoRisk"],    #list of 15
                 smooth = F)       # 曲线是否光滑，当光滑时，无法计算置信区间 
  # 计算AUC值
  auc<-auc(rocobj1)[1]
  auc_table=data.frame(Gene="nomoRisk",AUC=auc)
  auc_text<-paste0("AUC = ",round(auc,4))
  # 绘图
  p_roc=ggroc(rocobj1,
              color=colors[1],
              size=1,
              legacy.axes = F # FALSE时 横坐标为1-0 specificity；TRUE时 横坐标为0-1 1-specificity
  )+
    theme_classic()+
    geom_ribbon(aes(x=specificity, ymin = 0, ymax = sensitivity),
                fill = colors[1], alpha = .6)+
    # geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),        # 绘制对角线
    #              colour='grey', 
    #              linetype = 'dotdash'
    # ) +
    annotate('text',x = 0.25,y = 0.25,label=paste0("nomoRisk",'\n',auc_text),size=2.5)+
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  pdf(paste0("Nomogram/plot/ROC_valid_use.pdf"),width = 5.1653/2.54,height = 4.5991/2.54)
  print(p_roc)
  while (!is.null(dev.list()))  dev.off()
  pdf(paste0("Nomogram/plot/ROC_valid.pdf"),width = 6/2.54,height = 6/2.54)
  print(p_roc)
  while (!is.null(dev.list()))  dev.off()
}
