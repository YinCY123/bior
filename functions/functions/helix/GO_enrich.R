####GO富集####
GO_enrich<-function(genes=NULL,sp="org.Hs.eg.db",gene_type="SYMBOL",
                    sig=FALSE,DEG=NULL,
                    colors=c('#4DBBD5','#E64B35','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148','#B09C85')){
  if(is.null(genes))stop("缺少genes")
  dir.create("GO_KEGG")
  ####id转换####
  genes_id <- bitr(genes,fromType = gene_type,toType = "ENTREZID",OrgDb = sp,drop = T)
  
  ####总####
  ego_n <- enrichGO(gene = genes_id$ENTREZID,
                    OrgDb = sp,
                    keyType = 'ENTREZID',
                    ont = "ALL", 
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = T)
  ego_list <- data.frame(ego_n)
  
  ####输出数据####
  dir.create("GO_KEGG/data")
  GO_list=list(ego_list=ego_list,ego_n=ego_n,
               genes=genes,
               genes_id=genes_id,
               sp=sp)
  save(GO_list,file = "GO_KEGG/data/GO_list.rda")
  if(!nrow(ego_list)>0)stop("未富集到任何一条通路")
  dir.create("GO_KEGG/Table")
  write.csv(ego_list,"GO_KEGG/Table/enrichGO.csv")
  ####绘图####
  dir.create("GO_KEGG/plot")
  GO_plot_list<-list()
  g<-ego_list
  g$LOG10pvalue <- -log10(g$p.adjust)
  g1=g %>% group_by(ONTOLOGY) %>% slice_max(n = 3,order_by = LOG10pvalue,with_ties = F)
  p<-ggdotchart(g1, x="Description", y="LOG10pvalue", color = "ONTOLOGY",
             palette = colors,
             sorting = "descending",   #上升排序，区别于desc
             add = "segments",    #增加线段
             xlab = 'GO',ylab = "-log10(pvalue)",
             rotate = T,       #横向显示
             dot.size = 4.5,        #圆圈大小
             ggtheme = theme_pubr()+
               theme(plot.title = element_blank(),plot.subtitle = element_blank(),
                     plot.background = element_blank(),plot.margin = margin(t=2,r=10,b=2,l=2,unit="pt"),
                     #panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
                     axis.ticks = element_line(size = 0.5),axis.line = element_line(size = 0.5),
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     axis.text.y = element_text(size = 8,colour = "black"),axis.text.x = element_text(size = 6,colour = "black"),
                     axis.title.x = element_text(size = 7,colour = "black"),axis.title.y = element_blank(),
                     legend.title = element_blank(),legend.text = element_text(size = 7,colour = "black"),
                     legend.background = element_blank(),legend.key = element_blank(),legend.position = "top",
                     legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),legend.box.spacing = unit(3,"pt"))
  )+scale_x_discrete(labels=function(x) str_wrap(x, width=35))
  
  pdf("GO_KEGG/plot/GO_dotchart.pdf",width = 8/2.54,height = 10/2.54)
  print(p)
  dev.off()
  GO_plot_list[["GO_dotchart"]]=p
  
  p<-enrich_plot(dat=g,groups = NULL, group_name = "ONTOLOGY",enrich_geom="Nor",
                colors = colors,ay = "LOG10pvalue", ax = "Description", selCol = "pvalue", top_n = 5,
                xlab="Description", ylab="-log10(pvalue)", flab="",flip = T ,facet=T,
                base_size = 6, text_x_angle = 0, 
                text_x_hjust = 1.01, text_x_vjust = 1.01, 
                legend_position = "top",sort="asc",
                bar_width = 0.7)
  
  pdf("GO_KEGG/plot/GO_bar1.pdf",width = 8/2.54,height = 10/2.54)
  print(p)
  dev.off()
  GO_plot_list[["GO_bar1"]]=p
  
  p<-enrich_plot(dat=g,groups = NULL, group_name = "ONTOLOGY",enrich_geom="Sub",
                 colors = colors,ay = "LOG10pvalue", ax = "Description", selCol = "pvalue", top_n = 5,
                 xlab="Description", ylab="-log10(pvalue)", flab="",flip = T ,facet=T,
                 base_size = 6, text_x_angle = 0, main_cat_size = 7,sub_cat_size = 6,
                 text_x_hjust = 1.01, text_x_vjust = 1.01, 
                 legend_position = "top",sort="desc",
                 bar_width = 0.7)+theme(axis.title.x = element_text(size = 7,colour = "black"),
                                        axis.title.y = element_text(size = 7,colour = "black"),
                                        axis.text.x = element_text(size = 6,colour = "black"))
  
  pdf("GO_KEGG/plot/GO_bar2.pdf",width = 9/2.54,height = 13/2.54)
  print(p)
  dev.off()
  GO_plot_list[["GO_bar2"]]=p
  
  p<-barplot(ego_n)
  pdf("GO_KEGG/plot/GO_bar3.pdf",width = 8/1.5,height = 10/1.5)
  print(p)
  dev.off()
  GO_plot_list[["GO_bar3"]]=p
  ####上下调####
  if(sig){

    if(is.null(DEG))stop("选择了带有上下调信息的富集分析，但缺少DEG表格")
    if(!all(c("name","logFC","group")%in%colnames(DEG)))stop("DEG中应当包含name、logFC、group")
    
    genes_id_up<-genes_id[genes_id[,gene_type]%in%DEG$name[DEG$group=='up'],]
    genes_id_down<-genes_id[genes_id[,gene_type]%in%DEG$name[DEG$group=='down'],]
    
    if(length(genes_id_up[,1])+length(genes_id_down[,1])==0)stop("输入的基因集不是DEG中的差异基因")
    
    ego_up <- enrichGO(gene = genes_id_up$ENTREZID,
                       OrgDb = sp,
                       keyType = 'ENTREZID',
                       ont = "ALL", 
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = T)
    ego_up_list <- data.frame(ego_up)
    
    ego_down <- enrichGO(gene = genes_id_down$ENTREZID,
                         OrgDb = sp,
                         keyType = 'ENTREZID',
                         ont = "ALL", 
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = T)
    ego_down_list <- data.frame(ego_down)
    
    ####输出数据####
    GO_list[["genes_id_up"]]=genes_id_up
    GO_list[["genes_id_down"]]=genes_id_down
    GO_list[["ego_up_list"]]=ego_up_list
    GO_list[["ego_down_list"]]=ego_down_list
    GO_list[["ego_up"]]=ego_up
    GO_list[["ego_down"]]=ego_down
    save(GO_list,file = "GO_KEGG/data/GO_list.rda")
    write.csv(ego_up_list,"GO_KEGG/Table/enrichGO_up.csv")
    write.csv(ego_down_list,"GO_KEGG/Table/enrichGO_down.csv")
    ####绘图####
    #up####
    if(nrow(ego_up_list)>0){
      dir.create("GO_KEGG/plot/up")
      g<-ego_up_list
      g$LOG10pvalue <- -log10(g$p.adjust)
      g1=g %>% group_by(ONTOLOGY) %>% slice_max(n = 3,order_by = LOG10pvalue,with_ties = F)
      p<-ggdotchart(g1, x="Description", y="LOG10pvalue", color = "ONTOLOGY",
                    palette = colors,
                    sorting = "descending",   #上升排序，区别于desc
                    add = "segments",    #增加线段
                    xlab = 'GO',ylab = "-log10(pvalue)",
                    rotate = T,       #横向显示
                    dot.size = 4.5,        #圆圈大小
                    ggtheme = theme_pubr()+
                      theme(plot.title = element_blank(),plot.subtitle = element_blank(),
                            plot.background = element_blank(),plot.margin = margin(t=2,r=10,b=2,l=2,unit="pt"),
                            #panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
                            axis.ticks = element_line(size = 0.5),axis.line = element_line(size = 0.5),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text.y = element_text(size = 8,colour = "black"),axis.text.x = element_text(size = 6,colour = "black"),
                            axis.title.x = element_text(size = 7,colour = "black"),axis.title.y = element_blank(),
                            legend.title = element_blank(),legend.text = element_text(size = 7,colour = "black"),
                            legend.background = element_blank(),legend.key = element_blank(),legend.position = "top",
                            legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),legend.box.spacing = unit(3,"pt"))
      )+scale_x_discrete(labels=function(x) str_wrap(x, width=35))
      
      pdf("GO_KEGG/plot/up/GO_dotchart_up.pdf",width = 8/2.54,height = 10/2.54)
      print(p)
      dev.off()
      GO_plot_list[["GO_dotchart_up"]]=p
      
      p<-enrich_plot(dat=g,groups = NULL, group_name = "ONTOLOGY",enrich_geom="Nor",
                     colors = colors,ay = "LOG10pvalue", ax = "Description", selCol = "pvalue", top_n = 5,
                     xlab="Description", ylab="-log10(pvalue)", flab="",flip = T ,facet=T,
                     base_size = 6, text_x_angle = 0, 
                     text_x_hjust = 1.01, text_x_vjust = 1.01, 
                     legend_position = "top",sort="asc",
                     bar_width = 0.7)
      
      pdf("GO_KEGG/plot/up/GO_bar1_up.pdf",width = 8/2.54,height = 10/2.54)
      print(p)
      dev.off()
      GO_plot_list[["GO_bar1_up"]]=p
      
      p<-enrich_plot(dat=g,groups = NULL, group_name = "ONTOLOGY",enrich_geom="Sub",
                     colors = colors,ay = "LOG10pvalue", ax = "Description", selCol = "pvalue", top_n = 5,
                     xlab="Description", ylab="-log10(pvalue)", flab="",flip = T ,facet=T,
                     base_size = 6, text_x_angle = 0, main_cat_size = 7,sub_cat_size = 6,
                     text_x_hjust = 1.01, text_x_vjust = 1.01, 
                     legend_position = "top",sort="desc",
                     bar_width = 0.7)+theme(axis.title.x = element_text(size = 7,colour = "black"),
                                            axis.title.y = element_text(size = 7,colour = "black"),
                                            axis.text.x = element_text(size = 6,colour = "black"))
      
      pdf("GO_KEGG/plot/up/GO_bar2_up.pdf",width = 9/2.54,height = 13/2.54)
      print(p)
      dev.off()
      GO_plot_list[["GO_bar2_up"]]=p
      
      p<-barplot(ego_up)
      pdf("GO_KEGG/plot/up/GO_bar3_up.pdf",width = 8/1.5,height = 10/1.5)
      print(p)
      dev.off()
      GO_plot_list[["GO_bar3_up"]]=p
    }
    #down####
    if(nrow(ego_down_list)>0){
      dir.create("GO_KEGG/plot/down")
      g<-ego_down_list
      g$LOG10pvalue <- -log10(g$p.adjust)
      g1=g %>% group_by(ONTOLOGY) %>% slice_max(n = 3,order_by = LOG10pvalue,with_ties = F)
      p<-ggdotchart(g1, x="Description", y="LOG10pvalue", color = "ONTOLOGY",
                    palette = colors,
                    sorting = "descending",   #上升排序，区别于desc
                    add = "segments",    #增加线段
                    xlab = 'GO',ylab = "-log10(pvalue)",
                    rotate = T,       #横向显示
                    dot.size = 4.5,        #圆圈大小
                    ggtheme = theme_pubr()+
                      theme(plot.title = element_blank(),plot.subtitle = element_blank(),
                            plot.background = element_blank(),plot.margin = margin(t=2,r=10,b=2,l=2,unit="pt"),
                            #panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
                            axis.ticks = element_line(size = 0.5),axis.line = element_line(size = 0.5),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text.y = element_text(size = 8,colour = "black"),axis.text.x = element_text(size = 6,colour = "black"),
                            axis.title.x = element_text(size = 7,colour = "black"),axis.title.y = element_blank(),
                            legend.title = element_blank(),legend.text = element_text(size = 7,colour = "black"),
                            legend.background = element_blank(),legend.key = element_blank(),legend.position = "top",
                            legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),legend.box.spacing = unit(3,"pt"))
      )+scale_x_discrete(labels=function(x) str_wrap(x, width=35))
      
      pdf("GO_KEGG/plot/down/GO_dotchart_down.pdf",width = 8/2.54,height = 10/2.54)
      print(p)
      dev.off()
      GO_plot_list[["GO_dotchart_down"]]=p
      
      p<-enrich_plot(dat=g,groups = NULL, group_name = "ONTOLOGY",enrich_geom="Nor",
                     colors = colors,ay = "LOG10pvalue", ax = "Description", selCol = "pvalue", top_n = 5,
                     xlab="Description", ylab="-log10(pvalue)", flab="",flip = T ,facet=T,
                     base_size = 6, text_x_angle = 0, 
                     text_x_hjust = 1.01, text_x_vjust = 1.01, 
                     legend_position = "top",sort="asc",
                     bar_width = 0.7)
      
      pdf("GO_KEGG/plot/down/GO_bar1_down.pdf",width = 8/2.54,height = 10/2.54)
      print(p)
      dev.off()
      GO_plot_list[["GO_bar1_down"]]=p
      
      p<-enrich_plot(dat=g,groups = NULL, group_name = "ONTOLOGY",enrich_geom="Sub",
                     colors = colors,ay = "LOG10pvalue", ax = "Description", selCol = "pvalue", top_n = 5,
                     xlab="Description", ylab="-log10(pvalue)", flab="",flip = T ,facet=T,
                     base_size = 6, text_x_angle = 0, main_cat_size = 7,sub_cat_size = 6,
                     text_x_hjust = 1.01, text_x_vjust = 1.01, 
                     legend_position = "top",sort="desc",
                     bar_width = 0.7)+theme(axis.title.x = element_text(size = 7,colour = "black"),
                                            axis.title.y = element_text(size = 7,colour = "black"),
                                            axis.text.x = element_text(size = 6,colour = "black"))
      
      pdf("GO_KEGG/plot/down/GO_bar2_down.pdf",width = 9/2.54,height = 13/2.54)
      print(p)
      dev.off()
      GO_plot_list[["GO_bar2_down"]]=p
      
      p<-barplot(ego_down)
      pdf("GO_KEGG/plot/down/GO_bar3_down.pdf",width = 8/1.5,height = 10/1.5)
      print(p)
      dev.off()
      GO_plot_list[["GO_bar3_down"]]=p
    }
    #all####
    if(nrow(ego_list)>4){
      genes<-DEG[DEG$name%in%genes_id[,gene_type],]
      genes1 <- genes[,c('name','logFC')]
      colnames(genes1) <- c('ID','logFC')
      
      ego_list1 <- ego_list[,c('ONTOLOGY','ID','Description','Count','p.adjust','geneID')]
      colnames(ego_list1) <- c('category','ID','term','count','adj_pval','genes')
      ego_list1$genes <- gsub("/",",",ego_list1$genes)
      ##计算z-score
      circ <- circle_dat(ego_list1, genes1)
      write.csv(circ,"GO_KEGG/Table/GO_zscore.csv",row.names = FALSE)
      x <- circ[,c(2,8)]
      x <- x[!duplicated(x),]
      circ1 <- merge(ego_list1,x,by = "ID")
      #选取子集
      circ2 <- circ1 %>% arrange(adj_pval) %>%  group_by(category) %>% do(head(.,n = 3))
      
      #气泡图
      Bubble_labels<-unique(-log10(circ$adj_pval))%>%sort(decreasing = T)%>%head(.,n=4)
      p<-myGOBubble(circ, title = '', colour = colors[1:3], display = 'single', labels = Bubble_labels[4],table.legend = F)
      while (!is.null(dev.list()))  dev.off()
      pdf("GO_KEGG/plot/GO_Bubble1.pdf",width = 6,height = 5)
      print(p)  
      dev.off()
      GO_plot_list[["GO_Bubble1"]]=p
      #圈图
      while (!is.null(dev.list()))  dev.off()
      pdf("GO_KEGG/plot/GO_Circle.pdf",width = 6,height = 3.4)
      print(myGOCircle(circ, table.legend = T,nsub = circ2$ID,label.fontface = "plain",label.size = 2))
      dev.off()
      #弦图
      circ3<-circ[circ$ID%in%circ2$ID,]
      top_ngene=2
      top2 = circ3 %>% group_by(ID) %>% slice_max(n = top_ngene,order_by = abs(logFC),with_ties = F)
      while(length(unique(top2$genes))<10){
        top_ngene=top_ngene+1
        top2 = circ3 %>% group_by(ID) %>% slice_max(n = top_ngene,order_by = abs(logFC),with_ties = F)
      print(top_ngene)}
      genes3<-genes1[genes1$ID%in%top2$genes,]
      chord <- mychord_dat(data = circ3, genes = genes3, process = circ2$ID)
      p<-myGOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,
                lfc.min = min(chord[,"logFC"]),lfc.max = max(chord[,"logFC"]),
                limit = c(0, 0),ribbon.col=colorRampPalette(brewer.pal(12, "Set3"))(length(circ2$ID)))
      while (!is.null(dev.list()))  dev.off()
      pdf("GO_KEGG/plot/GO_Chord.pdf",width = 6,height = 7)
      print(p)
      dev.off()
      GO_plot_list[["GO_Chord"]]=p
    }else print("未富集到四个以上的通路，不绘制圈图、弦图、气泡图")
  }
  save(GO_plot_list,file="GO_KEGG/data/GO_plot_list.rda")
}




