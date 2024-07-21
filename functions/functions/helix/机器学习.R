
library(survival)
#install.packages("randomForestSRC")
library(randomForestSRC)
library(glmnet)
#install.packages("plsRcox")
library(plsRcox)
#install.packages("superpc")
library(superpc)
#install.packages("gbm")
library(gbm)
#devtools::install_github("binderh/CoxBoost")
library(CoxBoost)
#install.packages("survivalsvm")
library(survivalsvm)
library(dplyr)
library(tibble)
#install.packages("BART")
library(BART)
#### 参数检测 ####

if(is.null(ml_data))stop("缺少主要数据")
if(is.null(seed))stop("缺少随机数种子")
if(is.null(colors))stop("缺少颜色参数")
if(is.null(survminer_cut))stop("缺少阈值参数")
if(is.null(years))stop("缺少时间参数")

#调整颜色顺序
colors_red=(strtoi(stringr::str_sub(stringr::str_sub(colors,2),1,2),16)-strtoi(stringr::str_sub(stringr::str_sub(colors,2),5,6),16))
if(colors_red[1]>colors_red[2])colors[c(1,2)]=colors[c(2,1)]
#### 准备工作 ####
dir.create("Machine_learning")
ml_data <- lapply(ml_data,function(x){
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})

result <- data.frame()
est_data <- ml_data[[1]]
val_data_list <- ml_data
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',pre_var)]})

rf_nodesize <- 5#分类默认1，回归默认5

result_rs<-list()
result_gene<-list()
result_model<-list()
#### 1-1.RSF ####

set.seed(seed)
rf_select=var.select(Surv(OS.time,OS)~.,data = est_dd,
                     ntree = 1000,nodesize = rf_nodesize,
                     splitrule = 'logrank',
                     importance = T,
                     proximity = T,
                     forest = T,conservative = ("high"),
                     seed = seed)
if(length(rf_select$topvars)>5){
  est_dd_rf <- est_data[,c('OS.time','OS',rf_select$topvars)]
}else est_dd_rf =est_dd

fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd_rf,
             ntree = 1000,nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)

result_model[["RSF"]]=fit
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
result_rs[["RSF"]]=rs
result_gene[["RSF"]]=fit$xvar.names
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result,cc)

#### 2-1.Enet ####

x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time,est_dd$OS))

for (alpha in seq(0,1,0.1)) {
  alpha1=ifelse(alpha==0,"Ridge",
                ifelse(alpha==1,"LASSO",paste0('Enet','[alpha=',alpha,']')))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  result_model[[alpha1]]=fit
  result_rs[[alpha1]]=rs
  coefficient <- coef(fit, s="lambda.min")
  Active.Index <- which(as.numeric(coefficient) != 0)
  result_gene[[alpha1]]=rownames(coefficient)[Active.Index]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- alpha1
  result <- rbind(result,cc)
}

#### 3-1.StepCox ####

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  result_model[[paste0('StepCox','[',direction,']')]]=fit
  result_rs[[paste0('StepCox','[',direction,']')]]=rs
  result_gene[[paste0('StepCox','[',direction,']')]]=names(fit$coefficients)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}

#### 3-8.StepCox+gbm ####

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
  
  fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
  result_model[[paste0('StepCox','[',direction,']',' + survival-SVM')]]=fit
  result_rs[[paste0('StepCox','[',direction,']',' + survival-SVM')]]=rs
  result_gene[[paste0('StepCox','[',direction,']',' + survival-SVM')]]=fit$var.names
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + survival-SVM')
  result <- rbind(result,cc)
}

#### 4-1.CoxBoost ####

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
result_model[["CoxBoost"]]=fit
result_rs[["CoxBoost"]]=rs
result_gene[["CoxBoost"]]=fit$xnames
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)

#### 5.plsRcox####

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd[,pre_var],time=est_dd$OS.time,status=est_dd$OS),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd[,pre_var],time=est_dd$OS.time,event=est_dd$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
result_model[["plsRcox"]]=fit
result_rs[["plsRcox"]]=rs
result_gene[["plsRcox"]]=names(fit$dataX)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)

#### 6.superpc####

data <- list(x=t(est_dd[,-c(1,2)]),y=est_dd$OS.time,censoring.status=est_dd$OS,featurenames=colnames(est_dd)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
result_model[["SuperPC"]]=cv.fit
result_rs[["SuperPC"]]=rs
result_gene[["SuperPC"]]=names(fit$feature.scores)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result,cc)

#### 7.GBM ####

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)

rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
result_model[["GBM"]]=fit
result_rs[["GBM"]]<-rs
result_gene[["GBM"]]<-fit$var.names
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)

#### 8.survivalsvm ####

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd, gamma.mu = 1)

rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
result_model[["survival-SVM"]]=fit
result_rs[["survival-SVM"]]=rs
result_gene[["survival-SVM"]]=fit$var.names
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('survival-SVM')
result <- rbind(result,cc)

####组合####
result2 <- result
library(ggplot2)
library(ggsci)
library(tidyr)
range(result2$Cindex)
result2%>%filter(ID!=names(val_data_list)[1])%>%
  ggplot(aes(Cindex,reorder(Model,Cindex)))+
  geom_bar(width = 0.7,stat = 'summary',fun='mean',fill='orange2')+
  theme_classic()+
  labs(y=NULL)


dd2 <- pivot_wider(result2,names_from = 'ID',values_from = 'Cindex')%>%as.data.frame()
dd2[,-1] <- apply(dd2[,-1],2,as.numeric)
#提出训练集
dd3=dd2[,names(val_data_list)[1]]
names(dd3)=dd2$Model

model_list=c("LASSO","StepCox[both]","StepCox[backward]","CoxBoost","RSF")
model_list1=dd2$Model
#变量筛选的模型与其他模型
for(i in model_list){
  for(j in model_list1){
    if(i==j)next
    if((grepl("Enet",i)|grepl("LASSO",i)|grepl("Ridge",i))&(grepl("Enet",j)|grepl("LASSO",j)|grepl("Ridge",j)))next
    if(grepl("StepCox",i)&(grepl("StepCox",j)|grepl("SVM",j)))next
    if(grepl("[+]",paste0(i,j)))next
    weight1=dd3[i]/(dd3[i]+dd3[j])
    weight2=dd3[j]/(dd3[i]+dd3[j])
    rs <- lapply(names(val_dd_list),function(x){cbind(val_dd_list[[x]][,1:2],RS=result_rs[[i]][[x]]$RS*weight1+result_rs[[j]][[x]]$RS*weight2)})
    names(rs)=names(val_data_list)
    result_rs[[paste0(i," + ",j)]]=rs
    result_gene[[paste0(i," + ",j)]]=intersect(result_gene[[i]],result_gene[[j]])
  }
  model_list1=model_list1[model_list1!=i]
}
#变量筛选的模型与变量筛选的模型
for(i in model_list){
  for(j in model_list){
    if(i==j)next
    if((grepl("Enet",i)|grepl("LASSO",i)|grepl("Ridge",i))&(grepl("Enet",j)|grepl("LASSO",j)|grepl("Ridge",j)))next
    if(paste0(i,"+",j)%in%names(result_rs))next
    if(grepl("StepCox",i)&grepl("StepCox",j))next
    if(grepl("[+]",paste0(i,j)))next
    weight1=dd3[i]/(dd3[i]+dd3[j])
    weight2=dd3[j]/(dd3[i]+dd3[j])
    rs <- lapply(names(val_dd_list),function(x){cbind(val_dd_list[[x]][,1:2],RS=result_rs[[i]][[x]]$RS*weight1+result_rs[[j]][[x]]$RS*weight2)})
    names(rs)=names(val_data_list)
    result_rs[[paste0(i," + ",j)]]=rs
    result_gene[[paste0(i," + ",j)]]=intersect(result_gene[[i]],result_gene[[j]])
  }
}

result_101=data.frame()
for(i in 1:length(result_rs)){
  cc <- data.frame(Cindex=sapply(result_rs[[i]],function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- names(result_rs)[i]
  result_101 <- rbind(result_101,cc)
}

####Cindex热图####
setwd("./Machine_learning")
dir.create("plot")
#数据准备
result_101_met<-tidyr::spread(result_101,key = "ID",value = "Cindex")%>%
  tibble::column_to_rownames(var = "Model")%>%as.matrix()
#提出验证集
result_101_met<-result_101_met[,colnames(result_101_met)!=names(val_data_list)[1]]
result_101_met=result_101_met[order(rowMeans(result_101_met),decreasing = T),]
library(ComplexHeatmap)
#颜色
col_fun=circlize::colorRamp2(c(0.5,0.65,0.8), c("#3CC2B3", "white", "#ECBA5B"))
col_data<-colors[1:ncol(result_101_met)]
names(col_data)=colnames(result_101_met)
#注释条
top_annotation=HeatmapAnnotation(datasets=colnames(result_101_met),show_annotation_name = F,
                                 col = list(datasets=col_data),
                                 annotation_legend_param = list(datasets=list(labels_gp=gpar(fontsize=7),
                                                                              title_gp=gpar(fontsize=7))))

anno_pct = function(x) {
  
  max_x = max(x)
  text = paste0(sprintf("%.3f", x))
  cell_fun_pct = function(i) {
    pushViewport(viewport(xscale = c(0, max_x)))
    grid.roundrect(x = unit(1, "pt"), width = unit(x[i], "native"), 
                   height = unit(7, "pt"), 
                   just = "left", gp = gpar(fill = "#0000FF80", col = NA))
    grid.text(text[i], x = unit(x[i], "native")-unit(16, "pt"), just = "left",gp = gpar(fontsize=6))
    popViewport()
  }
  
  AnnotationFunction(
    cell_fun = cell_fun_pct,
    var_import = list(max_x, x, text), 
    which = "row",
    width = max_text_width(text)*1.25
  )
}
right_annotation=rowAnnotation("Mean Cindex"=anno_pct(rowMeans(result_101_met)),
                               annotation_name_gp= gpar(fontsize = 7),
                               annotation_name_rot=0,
                               annotation_name_side="top")
#添加数值
cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.3f", result_101_met[i, j]), x, y, gp = gpar(fontsize = 7))
}
#绘图
ha<-Heatmap(matrix = result_101_met,name = "Cindex",cluster_columns = F,cluster_rows = F,
          rect_gp = gpar(col = "white", lwd = 1),column_title_gp = gpar(fontsize=0),
          heatmap_legend_param = list(title_gp=gpar(fontsize = 7),labels_gp=gpar(fontsize = 7)),
          row_names_gp = gpar(fontsize=7),show_column_names = F,
          column_split =colnames(result_101_met),
          top_annotation = top_annotation,right_annotation = right_annotation,
          cell_fun = cell_fun,
          col = col_fun)

pdf("./plot/Cindex_heatmap.pdf",width = (6+(ncol(result_101_met)*2))/2.54,height = 42/2.54)
draw(ha)
dev.off()

####风险分数全家桶####
setwd("plot")
library(timeROC)
library(survminer)

best_model<-rownames(result_101_met)[1]
best_model_rs<-result_rs[[best_model]]
best_model_gene<-result_gene[[best_model]]
cutoff<-surv_cutpoint(best_model_rs[[names(val_data_list)[1]]],time = "OS.time",event = "OS",variables = "RS",minprop = 0.3)
cut=cutoff$cutpoint$cutpoint
best_model_cut=ifelse(survminer_cut,cut,median(best_model_rs[[names(val_data_list)[1]]][,"RS"]))

best_model_info=data.frame()
for(i in names(val_data_list)){
#时间依赖性ROC曲线
  dir.create(i)
  rt=best_model_rs[[i]]
  #ROC曲线
  ROC_rt=timeROC(T=rt$OS.time, delta=rt$OS,
                 marker=rt$RS, cause=1,
                 weighting='aalen',
                 times=c(1,3,5), ROC=TRUE)
  pdf(file=paste0("./",i,"/ROC_",i,".pdf"),width=5,height=5)
  plot(ROC_rt,time=years[1],col=colors[3],title=FALSE,lwd=2)
  plot(ROC_rt,time=years[2],col=colors[2],add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=years[3],col=colors[1],add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at ',years[1],' years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at ',years[2],' years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at ',years[3],' years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c(colors[3],colors[2],colors[1]),lwd=2,bty = 'n')
  dev.off()

#K-M生存曲线

  #读取输入文件
  rt=best_model_rs[[i]]
  rt$risk=ifelse(rt$RS>=best_model_cut,"High","Low")
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(OS.time, OS) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.05){
    pValue="p<0.05"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS.time, OS) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=3,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c(colors[2], colors[1]),
                     risk.table=FALSE,
                     risk.table.title="",
                     risk.table.height=.25,
                     fontsize = 2,
                     ggtheme = theme(axis.title = element_text(size = 7),
                                     axis.text = element_text(size = 6,colour = "black"),
                                     panel.background = element_blank(),
                                     axis.line = element_line(size = 0.5),
                                     legend.key.size = unit(0.5,"cm"),
                                     legend.text = element_text(size = 7),
                                     legend.title = element_blank())
  )
  while (!is.null(dev.list()))  dev.off()
  pdf(file=paste0("./",i,"/KM_",i,".pdf"),onefile = FALSE,width = 8/2.54,height =6/2.54)
  par(family="sans",ps="7")
  print(surPlot)
  dev.off()

  nHigh=as.numeric(sum(rt$risk=="High"))
  nLow=as.numeric(sum(rt$risk=="Low"))
  pValue=1-pchisq(diff$chisq,df=1)
#风险三联图

  rt=cbind(rt[,c("RS","risk")],val_data_list[[i]][rownames(rt),])   #读取输入文件
  rt$RS[rt$RS>quantile(rt$RS,0.99)]=quantile(rt$RS,0.99)
  rt$risk=factor(rt$risk, levels=c("Low", "High"))
  rt=rt[order(rt$RS),]      #按照风险得分对样品排序
  
  #绘制风险曲线
  lowMax=max(rt$RS[rt$risk=="Low"])
  line=rt[,"RS"]
  pdf(file=paste0("./",i,"/riskScore_",i,".pdf"), width=6, height=5)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep(colors[1],nLow),rep(colors[2],nHigh)))
  abline(h=lowMax,v=nLow,lty=2)
  legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c(colors[2],colors[1]),cex=1.2)
  dev.off()
  
  #绘制生存状态图
  color=as.vector(rt$OS)
  color[color==1]=colors[2]
  color[color==0]=colors[1]
  pdf(file=paste0("./",i,"/survStat_",i,".pdf"), width=6, height=5)
  plot(rt$OS.time, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c(colors[2],colors[1]),cex=1.2)
  abline(v=nLow,lty=2)
  dev.off()
  
  #绘制热图
  rt1 <- rt[,best_model_gene]
  mat=t(rt1)
  mat[mat > 3] = 3 #限定上限，使表达量大于3的等于3
  mat[mat < -3] = -3 #限定下限，使表达量小于-3的等于-3
  color <- colorRampPalette(c(colors[1],'white',colors[2]))(10)
  p <- pheatmap(mat,show_rownames = T,show_colnames = F,
                cluster_rows = F,cluster_cols = F,
                scale =F,fontsize = 7,name = " ",
                color = color)
  pdf(file=paste0("./",i,"/geneHeatmap_",i,".pdf"),width = length(best_model_gene)/4*1.4*1.4,height = length(best_model_gene)/4)
  print(p)
  dev.off()
  cor_info=data.frame(t(apply(mat,1,function(x)cor(x,rt$RS,method = "spearman"))))

  info=data.frame(dataset=i,
                  AUC1=sprintf("%.03f",ROC_rt$AUC[1]),
                  AUC2=sprintf("%.03f",ROC_rt$AUC[2]),
                  AUC3=sprintf("%.03f",ROC_rt$AUC[3]),
                  p=pValue,
                  nHigh=nHigh,
                  nLow=nLow)
  colnames(info)[2:4]=paste0("AUC_years",years)
  info=cbind(info,cor_info)
  best_model_info=rbind(best_model_info,info)
}
setwd("../")

####输出文件####
dir.create("data")
all_list<-list(val_data_list=val_data_list,
               base_model=result_model,
               all_gene=result_gene,
               all_Cindex=result_101,
               all_rs=result_rs,
               Cindex_heatmap_mat=result_101_met)
save(all_list,file = "./data/all_list.rda")

best_model_list<-list(best_model=best_model,
                      best_model_gene=best_model_gene,
                      best_model_cut=best_model_cut,
                      best_model_rs=best_model_rs,
                      best_model_info=best_model_info)
save(best_model_list,file = "./data/best_model_list.rda")

dir.create("Table")
write.csv(best_model_gene,"./Table/best_model_gene.csv")

setwd("../")
