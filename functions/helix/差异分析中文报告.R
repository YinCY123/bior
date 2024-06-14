load("DEG/data/DEG_list.rda")

####差异分析中文报告####

setwd("./DEG/")
dir.create('Report')
cn_path="/opt/scripts_repo/bin/template_cn.docx"

#报告文档
#officer
#word
#https://www.rdocumentation.org/packages/officer/versions/0.5.2
pkgs <- c('officer')
if(!requireNamespace("officer", quietly = TRUE))install.packages("officer")
library(officer)
library(magrittr)#只有加载了这个包才能用管道函数
library(tidyverse)

#创建空的word表格
# 
# my_doc<- read_docx(cn_path)

#了解常见排版的参数，这一步可以省略，只是为了更加清楚有哪些参数而已。
# s=styles_info(my_doc)
# test <- styles_info(my_doc)#比较重要的有：style_type,style_name
#设置字体
# text_font=fp_text(color = 'black',font.size = 12, font.family ='宋体') 
# title_font=fp_text(color = 'black',font.size = 15, font.family = '宋体',bold = T)
underline_font=fp_text(color = 'black',shading.color='yellow',
                       font.family = "Times New Roman",hansi.family = 'Times New Roman',
                       eastasia.family = "宋体",font.size = 10.5)

#写入内容
####材料与方法####
if(DEG_list$method=="edgeR"){
  Methods<-read_docx(cn_path) %>%        
    body_remove()%>%  
    body_add_par('材料与方法',style = "标题2") %>% 
    body_add_fpar(fpar(ftext('差异分析')),style = "标题3") %>% 
    body_add_fpar(fpar(ftext('本研究使用R包“edgeR（版本 3.36.0）”[PMID: 24753412]根据count矩阵鉴定对照组（'),
                       ftext(DEG_list$control,prop=underline_font),
                       ftext('）与实验组（'),
                       ftext(DEG_list$case,prop=underline_font),
                       ftext(paste0('）之间的差异',ifelse(DEG_list$feature_type=="基因",'表达基因（DEGs）',feature_type),'，筛选标准为')),
                       ftext(paste0("|log2Fold Change| >",DEG_list$logFC_value,'和',ifelse(DEG_list$P=="p.adj","FDR矫正后的P值 <","P值 <"),DEG_list$P_value,"。"),prop = underline_font)),
                  style = "Normal")%>% 
    print(target = "Report/Methods.docx")
}else if(DEG_list$method=="DESeq2"){
  Methods<-read_docx(cn_path) %>%        
    body_remove()%>%  
    body_add_par('材料与方法',style = "标题2") %>% 
    body_add_fpar(fpar(ftext('差异分析')),style = "标题3") %>% 
    body_add_fpar(fpar(ftext('本研究使用R包“DESeq2（版本 1.34.0）”[PMID: 25516281]根据count矩阵鉴定对照组（'),
                       ftext(DEG_list$control,prop=underline_font),
                       ftext('）与实验组（'),
                       ftext(DEG_list$case,prop=underline_font),
                       ftext(paste0('）之间的差异',ifelse(DEG_list$feature_type=="基因",'表达基因（DEGs）',feature_type),'，筛选标准为')),
                       ftext(paste0("|log2Fold Change| >",DEG_list$logFC_value,'和',ifelse(DEG_list$P=="p.adj","FDR矫正后的P值 <","P值 <"),DEG_list$P_value,"。"),prop = underline_font)),
                  style = "Normal")%>% 
    print(target = "Report/Methods.docx")
}else{
  Methods<-read_docx(cn_path) %>%        
    body_remove()%>%  
    body_add_par('材料与方法',style = "标题2") %>% 
    body_add_fpar(fpar(ftext('差异分析')),style = "标题3") %>% 
    body_add_fpar(fpar(ftext('本研究使用R包“limma （版本 3.50.0）”[PMID: 25605792]鉴定对照组（'),
                       ftext(DEG_list$control,prop=underline_font),
                       ftext('）与实验组（'),
                       ftext(DEG_list$case,prop=underline_font),
                       ftext(paste0('）之间的差异',ifelse(DEG_list$feature_type=="基因",'表达基因（DEGs）',feature_type),'，筛选标准为')),
                       ftext(paste0("|log2Fold Change| >",DEG_list$logFC_value,'和',ifelse(DEG_list$P=="p.adj","BH矫正后的P值 <","P值 <"),DEG_list$P_value,"。"),prop = underline_font)),
                  style = "Normal")%>% 
    print(target = "Report/Methods.docx")
}


####结果####

top10_padj = DEG_list$DEG[DEG_list$DEG$group!='no',] %>% group_by(group) %>% slice_max(n = 5,order_by = -log10(P.adj),with_ties = F)#选择P值最小的差异基因上下调各五个
top10_p = DEG_list$DEG[DEG_list$DEG$group!='no',] %>% group_by(group) %>% slice_max(n=5,order_by = -log10(P.value),with_ties = F)
if(DEG_list$P=="p")label=top10_p else label=top10_padj#选择P值最小的差异基因上下调各五个

Result<-read_docx(cn_path) %>%        
  body_remove()%>%  
  body_add_par('结果',style = "标题2") %>% 
  body_add_fpar(fpar(ftext(paste0("差异",DEG_list$feature_type,"分析"))),style = "标题3") %>% 
  body_add_fpar(fpar(ftext('通过对'),
                     ftext(DEG_list$case,prop=underline_font),
                     ftext('和'),
                     ftext(DEG_list$control,prop=underline_font),
                     ftext('的比较，共确定了'),
                     ftext(sum(DEG_list$DEG$group!="no"),prop=underline_font),
                     ftext(paste0("个差异",ifelse(DEG_list$feature_type=="基因",'表达基因（DEGs）',feature_type))),
                     ftext('，这些基因在两组之间的差异具有统计学意义（'),
                     ftext(paste0("|log2Fold Change| >",DEG_list$logFC_value,'，',ifelse(DEG_list$P=="p.adj","FDR矫正后的P值 <","P值 <"),DEG_list$P_value),prop = underline_font),
                     ftext('）。在'),
                     ftext(paste0(DEG_list$case,"样本中，",sum(DEG_list$DEG$group=="up"),"个基因上调，",sum(DEG_list$DEG$group=="down"),"个基因下调（Table）"),prop=underline_font),
                     ftext('。所有DEGs通过火山图'),
                     ftext("和热图",prop=underline_font),
                     ftext('进行可视化展示（Figure）。此外，还通过'),
                     ftext("热图和箱线图/小提琴图",prop=underline_font),
                     ftext('展示了P值排名前5位上调基因（'),
                     ftext(paste(label$name[label$group=="up"],collapse = "，"),prop=underline_font),
                     ftext('）和排名前5位下调基因（'),
                     ftext(paste(label$name[label$group=="down"],collapse = "，"),prop=underline_font),
                     ftext("）（Figure）。",prop = underline_font)),
                style = "Normal")%>%
  print(target = "Report/Result.docx")

####图注####
source("/opt/scripts_repo/bin/color_to_name.R")
color_name1=color_to_name(DEG_list$colors[1])
color_name2=color_to_name(DEG_list$colors[2])

Legend<-read_docx(cn_path) %>%        
  body_remove()%>%  
  body_add_par('图注合集',style = "标题2") %>% 
  body_add_fpar(fpar(ftext('火山图')),style = "标题3") %>% 
  body_add_fpar(fpar(ftext('火山图描述'),
                     ftext(DEG_list$case,prop=underline_font),
                     ftext('与'),
                     ftext(DEG_list$control,prop=underline_font),
                     ftext('样本之间差异基因的分布情况。'),
                     ftext(paste0(color_name2[2],"、",color_name1[2]),prop=underline_font),
                     ftext('和灰色圆点分别表示与上调、下调和无显著表达相关的基因表达水平。')),
                style = "图注")%>% 
  body_add_fpar(fpar(ftext('热图')),style = "标题3") %>% 
  body_add_fpar(fpar(ftext("（top10热图）热图描述了P值排名靠前的5个上调基因和5个下调基因的表达情况。")),
                style = "图注")%>% 
  body_add_fpar(fpar(ftext("（allgene热图）热图描述了所有差异基因的表达水平。")),
                style = "图注")%>% 
  body_add_fpar(fpar(ftext('箱线图')),style = "标题3") %>% 
  body_add_fpar(fpar(ftext('箱线图展示了P值排名靠前的5个上调基因和5个下调基因在'),
                     ftext(paste0(DEG_list$case,"与",DEG_list$control),prop=underline_font),
                     ftext('中的表达情况。')),
                style = "图注")%>% 
  body_add_fpar(fpar(ftext('小提琴图')),style = "标题3") %>% 
  body_add_fpar(fpar(ftext('小提琴图展示了P值排名靠前的5个上调基因和5个下调基因在'),
                     ftext(paste0(DEG_list$case,"与",DEG_list$control),prop=underline_font),
                     ftext('中的表达情况。')),
                style = "图注")%>% 
  print(target = "Report/legend.docx")

setwd("../")
