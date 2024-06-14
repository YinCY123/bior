dir.create("Cox")
dir.create("Cox/data")
dir.create("Cox/plot")
dir.create("Cox/Table")
cox_list=list(data=data,
              feature=feature,
              time_use=time_use,
              event_use=event_use,
              continuous=continuous,
              colors=colors,
              cox_type=cox_type)
####分析####
y<- Surv(time = data[,time_use],event = data[,event_use])#1为感兴趣事件
#批量单因素回归模型建立：Uni_cox_model Multi_cox_model
Uni_cox_model<-
  function(x){
    FML <- as.formula(paste0 ("y~",x))
    cox<- coxph(FML,data=data)
    cox1<-summary(cox)
    HR <- round(cox1$coefficients[,2],2)#提取HR值，保留2位小数
    PValue <- ifelse(cox1$coefficients[,5] < 0.05,"<0.05",round(cox1$coefficients[,5],3))#提取p值，保留3位小数
    CI5 <-round(cox1$conf.int[,3],2)#提取CI，保留2位小数
    CI95 <-round(cox1$conf.int[,4],2)
    #将提取到的信息放入表格中（Uni_cox_model）
    Uni_cox_model<- data.frame(
      names <-rownames(cox1$conf.int),#第1列为亚变量名
      'HR' = HR,#第2列为HR值
      'CI5' = CI5,#第3列为95%ci下区间
      'CI95' = CI95,#第4列为95%ci上区间
      'P' = PValue)#第5列为P值
    return(Uni_cox_model)#返回，开始，进行循环
  }  
#查看原始数据变量的名字
names(data)

if(cox_type=="Uni_cox"){
  Uni_cox <- lapply(feature, Uni_cox_model)
  library(plyr)
  Uni_cox <- ldply(Uni_cox,data.frame)
  #将95%CI连接起来
  Uni_cox$HR.CI95 <- paste0(Uni_cox$HR," (",Uni_cox$CI5,'-',Uni_cox$CI95,")");Uni_cox
  #第一列列名为'Characteristics'
  colnames(Uni_cox)[1] <- 'Characteristics'
  result <- Uni_cox[,c(1:4,6,5)]
}else if(cox_type=="Multi_cox"){
  FML=as.formula(paste0("y~",paste(feature,collapse = "+")))
  cox<- coxph(FML,data=data)
  cox1<-summary(cox)
  HR <- round(cox1$coefficients[,2],2)#提取HR值，保留2位小数
  PValue <- ifelse(cox1$coefficients[,5] < 0.05,"<0.05",round(cox1$coefficients[,5],3))#提取p值，保留3位小数
  CI5 <-round(cox1$conf.int[,3],2)#提取CI，保留2位小数
  CI95 <-round(cox1$conf.int[,4],2)
  #将提取到的信息放入表格中
  result<- data.frame(
    names <-rownames(cox1$conf.int),#第1列为亚变量名
    'HR' = HR,#第2列为HR值
    'CI5' = CI5,#第3列为95%ci下区间
    'CI95' = CI95,#第4列为95%ci上区间
    'P' = PValue)#第5列为P值
  result$HR.CI95 <- paste0(result$HR," (",result$CI5,'-',result$CI95,")")
  #第一列列名为'Characteristics'
  colnames(result)[1] <- 'Characteristics'
  result <- result[,c(1:4,6,5)]
}else stop("未输入正确的cox_type，请检查")

####表格整理
#删除部分变量名，只保留亚变量
result$Characteristics<-str_remove(result$Characteristics,"age|gender|stage|Tstage|Nstage|Mstage|risk")
result
cox_list[["result"]]=result
write.csv(result,paste0("Cox/Table/",cox_type,"_result.csv"))
#给参考变量插入空行
ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
#插入空行，形成一个新表
result1<-data.frame("Characteristics", NA, NA, NA, "HR(95%CI)","p")
colnames(result1)<-colnames(result)

if(!continuous){
  for(i in feature){
    a<-as.data.frame(table(data[,i]))
    result1<-rbind(result1,
                   ins(i),
                   ins(as.character(a$Var1[!a$Var1%in%result$Characteristics])),
                   result[which(feature==i),])
  }
}else{
  result1=rbind(result1,result)
}

result1<-rbind(result1,
               c(NA, NA, NA, NA, NA,NA))

rownames(result1)<-1:nrow(result1)

myVars <- feature
catVars <-  feature
library(tableone)
if(!continuous){
  table1<- print(CreateTableOne(vars=myVars,
                              data = data,
                              factorVars = catVars),
               showAllLevels=TRUE)

N<-data.frame(c(NA,NA),
              c(NA,NA))
colnames(N)<-colnames(table1)
for(i in 1:length(feature)){
  N=rbind(N,
          table1[(i*2):(i*2+1),],
          c(NA,NA))
  }  
  N<-N[,-1]
  N<-data.frame(N)
}else{
  N=data.frame(rep(NA,nrow(result1)))
  }


result2<-cbind(result1,N)
result2<-result2[,c(1,7,2:6)]

result2[1,]<-c("Characteristics","Number (%)",NA,NA,NA,"HR (95%CI)","P.value")

####图形进行美化
library(forestplot)
hrzl_lines<-list(gpar(lty=1,lwd=2),#表头上方添加实线
                 gpar(lty=2),#表头下方添加虚线
                 gpar(lwd=2,lty=1,columns=c(1:4)))#最后一行下方添加实线

names(hrzl_lines)<-c('1','2',nrow(result2)+1)

if(continuous){
  is_summary<-is.na(result2[,3]) #按顺序指定每行是否加粗，T为加粗，F不加粗
  is_summary[1]=TRUE
}else{
  is_summary<-is.na(result2[,2]) #按顺序指定每行是否加粗，T为加粗，F不加粗
  is_summary[1]=TRUE
  }

if(continuous)is_numeric=c(1,6,7) else is_numeric=c(1,2,6,7)

if(length(seq(min(result$CI5)-min(result$CI5)%%0.5,max(result$CI95),0.5))<=7){
  xticks=seq(min(result$CI5)-min(result$CI5)%%0.5,max(result$CI95)+0.5,0.5)
}else if(length(seq(min(result$CI5)-min(result$CI5)%%0.5,max(result$CI95),1))<=7){
  xticks=seq(min(result$CI5)-min(result$CI5)%%0.5,max(result$CI95)+1,1)
}else{
  xticks=seq(min(result$CI5)-min(result$CI5)%%0.5,max(result$CI95)+2,2)
}

fig2<- forestplot(result2[,is_numeric], 
                  mean=as.numeric(result2[,3]),   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                  lower=as.numeric(result2[,4]),  #告诉函数表格第3列为5%CI，
                  upper=as.numeric(result2[,5]),  #表格第5列为95%CI，它俩要化作线段，穿过方块
                  zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                  boxsize=0.25,       #设置小黑块的大小
                  graph.pos= "right" ,
                  hrzl_lines=hrzl_lines, #最后一行下方添加实线
                  graphwidth = unit(.25,"npc"),
                  xticks=xticks, #森林图刻度
                  is.summary=is_summary, #按顺序指定每行是否加粗，T为加粗，F不加粗
                  txt_gp=fpTxtGp(
                    label=gpar(cex=0.5),
                    ticks=gpar(cex=0.5), 
                    xlab=gpar(cex=0.5), 
                    title=gpar(cex=0.5)),
                  lwd.zero=0.5,
                  lwd.ci=1,
                  lwd.xaxis=2, 
                  lty.ci=1,
                  ci.vertices =T,
                  ci.vertices.height=0.2, 
                  clip=c(0,8),
                  ineheight=unit(9, 'mm'), 
                  line.margin=unit(9, 'mm'),
                  colgap=unit(1.7, 'mm'),
                  fn.ci_norm="fpDrawDiamondCI", 
                  title="",
                  col=fpColors(box =colors, 
                               lines =colors, 
                               zero = "black"))       #森林图应插在图形第2列
while (!is.null(dev.list()))  dev.off()
pdf(paste0("Cox/plot/",cox_type,".pdf"),height = 8/2.54,width = 8.5/2.54)
print(fig2)
dev.off()
cox_list[["xticks"]]=xticks
cox_list[["forestplot"]]=fig2
save(cox_list,file = paste0("Cox/data/",cox_type,"_list.rda"))