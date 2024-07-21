load("Machine_learning/data/best_model_list.rda")

####机器学习中文报告####

setwd("./Machine_learning/")
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

Methods<-read_docx(cn_path) %>%        
  body_remove()%>%  
  body_add_par('材料与方法',style = "标题2") %>% 
  body_add_fpar(fpar(ftext('基于机器学习的集成方法构建预后模型')),style = "标题3") %>% 
  body_add_fpar(fpar(ftext('为了开发具有高精度和稳定性能的'),
                     ftext(pheno_name,prop=underline_font),
                     ftext('相关预后模型，我们集成了10种机器学习算法，包括随机生存森林（RSF）、弹性网络（Enet）、Lasso、Ridge、StepCox、CoxBoost、Cox偏最小二乘回归（plsRcox）、监督主成分（SuperPC）、广义增强回归建模（GBM）和生存支持向量机（survival-SVM）。'),
                     ftext('一些算法具有特征选择的能力，例如 Lasso、StepCox、CoxBoost 和 RSF。因此，我们将这些算法结合起来生成了一个共识模型，总共进行了 101 种算法组合。'),
                     ftext('RSF模型是通过R包“randomForestSRC”实现的；Enet、Lasso和Ridge是通过R包“glmnet”实现的，L1-L2 权衡参数α设置为 0-1（间隔为0.1）；StepCox是通过R包“survival”实现的，逐步搜索的方向设为"both"、"backward"和"forward"；'),
                     ftext('CoxBoost是通过R包“CoxBoost”实现的，通过“optimCoxBoostPenalty”函数计算最佳惩罚参数，通过“cv.CoxBoost”函数计算最佳步骤数；plsRcox是通过R包“plsRcox”实现的，通过函数“cv.plsRcox”确定需要提取的组件数目；'),
                     ftext('SuperPC是通过R包“superpc”实现的，通过“superpc.cv”函数确定最佳阈值；GBM是通过R包“gbm”实现的，通过函数“cv.gbm”选择交叉验证误差最小的树的数量；survival-SVM是通过R包“survivalsvm”实现的。'),
                     ftext(paste0("模型的训练是在",names(best_model_list$best_model_rs)[1],'中进行的，并且在',paste(names(best_model_list$best_model_rs)[-1],collapse = "、")),prop = underline_font),
                     ftext('中进行验证。C-index值作为模型效能的评价指标，我们选择平均C-index值最高的模型作为最终结果。')),
                style = "Normal")%>% 
  print(target = "Report/Methods.docx")

####结果####

#AUC
for(i in 1:nrow(best_model_list$best_model_info)){
  if(i==1){
    temp_text=paste0(best_model_list$best_model_info[i,1],"数据集在",years[1],"年、",years[2],"年和",years[3],"年的 AUC 值分别为",best_model_list$best_model_info[i,2],"、",best_model_list$best_model_info[i,3],"和",best_model_list$best_model_info[i,4],"（Figure）")
    AUC_text=temp_text
  }else{
    temp_text=paste0(best_model_list$best_model_info[i,1],"数据集的 AUC 值分别为",best_model_list$best_model_info[i,2],"、",best_model_list$best_model_info[i,3],"和",best_model_list$best_model_info[i,4],"（Figure）")
    AUC_text=paste(AUC_text,temp_text,sep = "；")
  }
}

Result<-read_docx(cn_path) %>%        
  body_remove()%>%  
  body_add_par('结果',style = "标题2") %>% 
  body_add_fpar(fpar(ftext(paste0(pheno_name,"相关预后模型的构建"))),style = "标题3") %>% 
  body_add_fpar(fpar(ftext('我们采用了十种机器学习算法来构建'),
                     ftext(pheno_name,prop=underline_font),
                     ftext('相关预后模型。这些算法被应用于'),
                     ftext(paste0(names(best_model_list$best_model_rs)[1],'队列和',length(best_model_list$best_model_rs)-1,"个外部验证数据集（",paste(names(best_model_list$best_model_rs)[-1],collapse = "、"),"）","以确定最佳模型（",
                                  length(best_model_list$best_model_rs)-1,"个验证队列中平均C-index值最大的模型）。"),prop = underline_font),
                     ftext(paste0("最终的",best_model_list$best_model,"算法鉴定了",length(best_model_list$best_model_gene),"个最有价值的",pheno_name,"相关特征基因")),
                     ftext("（Table）",prop = underline_font),
                     ftext('构建了性能最佳的模型'),
                     ftext("（Figure）。",prop = underline_font),
                     ftext('值得注意的是，在'),
                     ftext(names(best_model_list$best_model_rs)[1],prop = underline_font),
                     ftext(paste0('中，随着风险分数的增加，',sum(best_model_list$best_model_info[1,best_model_list$best_model_gene]>0),'个特征基因表达升高，',sum(best_model_list$best_model_info[1,best_model_list$best_model_gene]<0),'个特征基因表达降低，并且死亡患者占比增加')),
                     ftext("（Figure S）。",prop = underline_font)),
                style = "Normal")%>%
  body_add_fpar(fpar(ftext('我们对'),
                     ftext(paste0(paste(names(best_model_list$best_model_rs)[-length(best_model_list$best_model_rs)],collapse = "、"),"和",names(best_model_list$best_model_rs)[length(best_model_list$best_model_rs)]),prop = underline_font),
                     ftext('数据集进行生存分析，高风险评分与生存时间缩短相关'),
                     ftext("（Figure）。",prop = underline_font),
                     ftext(AUC_text,prop = underline_font),
                     ftext(paste0("。这些结果强调了",pheno_name,"相关模型的预后意义。"))),
                style = "Normal")%>%
  print(target = "Report/Result.docx")

####图注####
#生存曲线
for(i in 1:nrow(best_model_list$best_model_info)){
  if(i==1){
    temp_text=paste0("KM生存曲线展示了",best_model_list$best_model_info[i,1],"数据集中风险评分与总生产期之间的联系")
    KM_text=temp_text
  }else{
    temp_text=paste0("KM生存曲线展示了",best_model_list$best_model_info[i,1],"数据集中风险评分与总生产期之间的联系")
    KM_text=paste(KM_text,temp_text,sep = "；")
  }
}

for(i in 1:nrow(best_model_list$best_model_info)){
  if(i==1){
    temp_text=paste0("ROC曲线展示了模型在",best_model_list$best_model_info[i,1],"中",years[1],"、",years[2],"、",years[3],"年的预后效能")
    ROC_text=temp_text
  }else{
    temp_text=paste0("ROC曲线展示了模型在",best_model_list$best_model_info[i,1],"中",years[1],"、",years[2],"、",years[3],"年的预后效能")
    ROC_text=paste(ROC_text,temp_text,sep = "；")
  }
}

Legend<-read_docx(cn_path) %>%        
  body_remove()%>%  
  body_add_par('图注合集',style = "标题2") %>% 
  body_add_fpar(fpar(ftext('Cindex热图')),style = "标题3") %>% 
  body_add_fpar(fpar(ftext('101种机器学习算法组合在不同数据集种的C-index值热图，右侧柱状图展示了算法组合在验证队列中的平均C-index值')),
                style = "图注")%>% 
  body_add_fpar(fpar(ftext('生存曲线')),style = "标题3") %>% 
  body_add_fpar(fpar(ftext(KM_text)),
                style = "图注")%>% 
  body_add_fpar(fpar(ftext('ROC曲线')),style = "标题3") %>% 
  body_add_fpar(fpar(ftext(ROC_text)),
                style = "图注")%>% 
  print(target = "Report/legend.docx")

setwd("../")
