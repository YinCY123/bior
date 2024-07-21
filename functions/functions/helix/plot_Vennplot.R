####Venn图####
library(grDevices)
library(ggforce)
plot_Vennplot=function(data=NULL,feature=NULL,colors=colors,feature_type=NULL,default_colors=TRUE){
  #数据检测####
  if(is.null(data))stop("缺少data")
  if(!is.null(feature)&!all(feature%in%names(data)))stop("feature不在data之中")
  if(is.null(feature_type))stop("缺少feature_type")
  #数据整理####
  if(!is.null(feature))data=data[feature]
  if(length(data)==1)stop("只输入了一个基因集")
  if(length(data)>length(colors))warning("colors数目不够，使用默认配色")
  dup_data=unlist(data)
  if(sum(duplicated(dup_data))==0)stop("所有集之间都没有交集")
  dup_data=as.data.frame(table(dup_data))
  all=sum(dup_data$Freq==length(data))
  if(all==0){
    Vennplot_text="没有共同交集，请自行提取关注的交集或者更换数据集"
    }else{
    all_data=as.character(dup_data$dup_data[dup_data$Freq==length(data)])
    all_text=ifelse(all>10,"（Table）",paste0("包括：",paste(all_data,collapse = "、")))
    Vennplot_text<-paste0(length(data),"个",feature_type,"集之间有",all,"个共同的交集",feature_type,all_text,"。")
  }
  Vennplot_list<<-list(Vennplot_text=Vennplot_text,Vennplot_intersect=all_data)
  #绘图####
  if(length(data)<=4){
    p<-ggvenn(data,
              fill_color = colors,fill_alpha = 0.4,
              show_percentage = F,set_name_size = 5,
              stroke_color = "white",stroke_size = 0.5)
  }else{
    if(all==0)stop("没有共同交集，请自行提取关注的交集或者更换数据集")
    data_ellipse=data.frame()
    n=length(data)
    a=5
    for(i in 1:n){
      temp=data.frame(x=cos((i-1)*(2*pi/n)+pi/2)*a,
                      y=sin((i-1)*(2*pi/n)+pi/2)*a,
                      angle=(i-1)*(2*pi/n)+pi/2)
      data_ellipse=rbind(data_ellipse,temp)
    }
    data_ellipse$name=names(data)
    data_ellipse$number=unlist(lapply(data,function(x)length(x)))
    data_ellipse$name_y=ifelse(abs(data_ellipse$y)<=1,-1,data_ellipse$y*2.5)
    p=ggplot(data=data_ellipse,aes(x0=x,y0=y,angle=angle,a=5,b=1.5,fill=name))+
      geom_ellipse(alpha=0.4)+
      geom_text(data = data_ellipse,aes(x=x*2.5,y=name_y,label=name,angle=0))+
      geom_text(data = data_ellipse,aes(x=x*1.5,y=y*1.5,label=number,angle=0))+
      scale_x_continuous(limits = c(-(5*2)-5,(5*2)+5))+
      geom_circle(aes(x0=0,y0=0,r=2),inherit.aes = F,fill="white")+
      annotate("text",x=0,y=0,label=as.character(all))+
      theme_classic()+
      theme(axis.ticks=element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),legend.position = "none")
    if(!default_colors){
      if(length(data)<=length(colors))p=p+scale_fill_manual(values = colors[1:length(data)])
    }
  }
  return(p)
}
