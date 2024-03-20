####数据整理####

rm(list = ls())#清空环境中所有变量
##导入文件
library(readxl)
AROW_CLINICAL <- read_excel("ROW-CLINICAL.xlsx", 
                            col_types = c("text", "numeric", "numeric", 
                                          "text", "text", "text", "text", "text", 
                                          "text", "text", "text", "text", "text", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric"))

# 改OS为int格式
AROW_CLINICAL$OS <- as.integer(AROW_CLINICAL$OS)
AROW_CLINICAL <- as.data.frame(AROW_CLINICAL)
AROW_CLINICAL <- AROW_CLINICAL[,-(4:5)]         ##去除RFS数据，只保存生存数据

#调整所有临床因素为数字形
rt <- AROW_CLINICAL
#设置参数     #从前往后分为为1,2，3,4,5    # 这里"female"会被编码为1，"male"会被编码为2
rt$GENDER <- as.numeric(factor(rt$GENDER,labels=c("Female", "Male")))
rt$T <- as.numeric(factor(rt$T,labels=c("T1", "T2", "T3", "T4", "TX")))
rt$N <- as.numeric(factor(rt$N,labels=c("N0", "N1", "NX")))
rt$M <- as.numeric(factor(rt$M,labels=c("M0", "M1", "MX")))
rt$HBV <- as.numeric(factor(rt$HBV,labels=c("True", "False")))
rt$HCV <- as.numeric(factor(rt$HCV,labels=c("True", "False")))
rt$GRADE <- as.numeric(factor(rt$GRADE,labels=c("G1", "G2", "G3", "G4", "GX")))
rt$ECGA <- as.numeric(factor(rt$ECGA,labels=c("E0", "E1", "E2", "E3", "E4", "EX")))
AROW_CLINICAL <- rt
rm(rt)

AROW_CLINICAL <- as.data.frame(AROW_CLINICAL)


#提取TCGA数据以及XIJING数据
ATCGA <- AROW_CLINICAL[-(1:50),]
AXIJING <- AROW_CLINICAL[(1:50),]
TCGA <- ATCGA

#####区分训练集和验证集####
#在TCGA中区分训练集和内部验证集，设置XIJING为外部验证集
library(caret)
set.seed(1241)    # 设置随机种子，以便结果可以复现
# createDataPartition返回的是行索引,,7:3
index <- createDataPartition(TCGA$OS, p = .7, list = FALSE)
# 根据索引创建训练集和测试集
ATCGA_TRAIN <- TCGA[index, ]
ATCGA_TEST <- TCGA[-index, ]
rm(TCGA)
rm(index)

#####  加入西京数据####
#套用别人代码，GSE57303 = ATCGA_TEST, GSE62254 = XIJING
mm <- list(TCGA = ATCGA_TRAIN,
           GSE57303 = ATCGA_TEST, GSE62254 = XIJING)

# 数据标准化
mm <- lapply(mm,function(x){
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})

result <- data.frame()
# TCGA作为训练集
est_data <- mm$TCGA
# GEO作为验证集
val_data_list <- mm

pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[, c('OS.time', 'OS', pre_var)]
val_dd_list <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', pre_var)]})





####加载包####
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(limma)
library(tidyverse)
library(dplyr)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

####设置随机数种子####
seed <- 12371

####随机森林模型构建####
#首先，定义一个评估模型性能的函数，这里以C指数为例####
evaluate_model <- function(est_data, val_dd_list, ntree, nodesize) {
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_data,
               ntree = ntree, nodesize = nodesize,
               splitrule = 'logrank',
               importance = TRUE, proximity = TRUE, forest = TRUE, seed = seed)
  
  rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
  cindex_values <- sapply(rs, function(x){
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
  })
  
  # 返回平均C指数作为性能评估
  mean_cindex <- mean(cindex_values)
  return(mean_cindex)
}

#接下来，进行网格搜索以寻找最佳的nodesize和ntree####
nodesize_options <- c(1, 5, 10, 15)  # 调整为需要测试的值
ntree_options <- c(500, 1000, 1500)  # 调整为需要测试的值

best_cindex <- 0
best_nodesize <- NA
best_ntree <- NA

for (nodesize in nodesize_options) {
  for (ntree in ntree_options) {
    current_cindex <- evaluate_model(est_dd, val_dd_list, ntree, nodesize)
    
    if (current_cindex > best_cindex) {
      best_cindex <- current_cindex
      best_nodesize <- nodesize
      best_ntree <- ntree
    }
  }
}

cat("最佳nodesize:", best_nodesize, "\n")
cat("最佳ntree:", best_ntree, "\n")
cat("最高C指数:", best_cindex, "\n")


####绘制最佳ntree图####
# 评估模型性能的函数，这里以C指数为例
evaluate_model_cindex <- function(est_data, val_dd_list, ntree) {
  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_data,
               ntree = ntree, seed = 130)
  
  rs <- lapply(val_dd_list, function(x){
    cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
  })
  cindex_values <- sapply(rs, function(x){
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
  })
  
  mean_cindex <- mean(cindex_values)
  return(mean_cindex)
}

# 定义ntree的范围
ntree_values <- seq(100, 1000, by = 100)
cindex_results <- numeric(length(ntree_values))

# 对每个ntree值评估模型性能
for (i in seq_along(ntree_values)) {
  cindex_results[i] <- evaluate_model_cindex(est_dd, val_dd_list, ntree_values[i])
}

# 绘图
plot_data <- data.frame(ntree = ntree_values, cindex = cindex_results)
ggplot(plot_data, aes(x = ntree, y = cindex)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Random Forest Model Performance vs. Number of Trees",
       x = "Number of Trees (ntree)",
       y = "C-index")



####重新进行随机森林处理####
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd,
                   ntree = best_ntree, nodesize = best_nodesize,
                   splitrule = 'logrank', importance = TRUE,
                   proximity = TRUE, forest = TRUE, seed = seed)


rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result, cc)


head(result)


####重要性导出并绘图####
# 假设 fit 是上述步骤中训练好的随机森林模型
# 导出变量重要性
var_importance <- fit[["importance"]]

# 将重要性数据转换为数据框，以便于绘图
importance_df <- data.frame(Variable = names(var_importance), Importance = var_importance)

# 按重要性排序
importance_df <- importance_df[order(-importance_df$Importance), ]

# 绘制重要性排序的柱状图
ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  coord_flip() +  # 将坐标轴翻转，使得变量名更容易阅读
  labs(title = "Variable Importance in Random Forest Model", x = "Variable", y = "Importance") +
  theme(plot.title = element_text(hjust = 0.5)) # 居中标题


####分别计算训练集的风险评分和c-index####






