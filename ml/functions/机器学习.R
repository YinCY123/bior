library(fs)
dir_create("20240528/res/fig5/Machine_learning")
# dir.create('Machine_learning')
exp_combat <- readRDS("20240528/res/mtx_combat.rds")
group_list<-readRDS("20240528/res/sample_info.rds")
if(!all(colnames(exp_combat)==group_list$title))stop('样本名未对应上')
hubgene <- read.table('20240528/res/fig3/WGCNA/ WGCNA RA glutamine_metabolism gene1.txt', header = F)[[1]]
input <- cbind(data.frame(group=ifelse(group_list$group=="normal",0,1)),t(exp_combat[hubgene,]))

timestart<-Sys.time()

####svm####

set.seed(101)
library(e1071)
source("20240528/functions/msvmRFE.R")

nfold = 10 #10倍交叉验证
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))

results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)
top.features = WriteFeatures(results, input, save=F)
featsweep = lapply(1:(ncol(input)-1), FeatSweep.wrap, results, input)

no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

pdf("20240528/res/fig5/Machine_learning/svm_rfe.pdf", height = 8, width = 10)
PlotErrors(errors, no.info=no.info)
dev.off()


#结束计时并输出运行时间
timeend<-Sys.time()
print(timeend-timestart)

#机器学习,svm-RFE。并行代码。Rmpi这个包遇到问题，无法多线程同时进行
# library(e1071)
# library(Rmpi)
# library(snow)
# library(parallel)
# 
# source("prepare/msvmRFE.R")
# 
# nfold = 10
# nrows = nrow(input)
# folds = rep(1:nfold, len=nrows)[sample(nrows)]
# folds = lapply(1:nfold, function(x) which(folds == x))
# 
# #make a cluster
# cl <- makeMPIcluster(mpi.universe.size())
# 
# clusterExport(cl, list("input","svmRFE","getWeights","svm"))
# results <-parLapply(cl,folds, svmRFE.wrap, input, k=10, halve.above=100)
# top.features = WriteFeatures(results, input, save=F)
# 
# clusterExport(cl, list("top.features","results", "tune","tune.control"))
# featsweep = parLapply(cl,1:100, FeatSweep.wrap, results, input)
# stopCluster(cl)
# 
# no.info = min(prop.table(table(input[,1])))
# errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
# #dev.new(width=4, height=4, bg='red')
# pdf("svm_rfe.pdf", height = 8, width = 10)
# PlotErrors(errors, no.info=no.info)
# dev.off()
# plot(top.features)
# mpi.exit()

####lasso####

library(glmnet)
# data("BinomialExample")
# x <- BinomialExample$x
# y <- BinomialExample$y
x <- input[,-1]
groups<-input$group
levels(groups) <- c(0, 1)
y <- as.integer(as.character(groups))
fit <- glmnet(x, y, family = "binomial", nlambda = 100, alpha = 1)
pdf("20240528/res/fig5/Machine_learning/lasso1.pdf", height = 8, width = 10)
plot(fit, xvar = "lambda", label = F)
dev.off()
cvfit <- cv.glmnet(data.matrix(x), y, nfolds = 10)
cvfit$lambda.min
cvfit$lambda.1se

pdf("20240528/res/fig5/Machine_learning/lasso2.pdf", height = 8, width = 10)
plot(cvfit)
dev.off()
coef <- coef(cvfit, s = "lambda.min")
factors <- as.matrix(coef)
factors[factors == 0] <- NA
factors <- na.omit(factors)
write.csv(factors, "20240528/res/fig5/Machine_learning/lasso.csv")

####randomforest####
#随机森林，行为样本，列为基因，搞一列为分类 

library(randomForest)

infr <- as.data.frame(input)
colnames(infr) <- make.names(colnames(infr))
infr[,1] <- as.factor(infr[,1])
model <- randomForest(group~., data = infr, proximity = TRUE, ntree = 1000)

#画图看种多少树合适
oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=3),
  Type=rep(c("OOB", control, case), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
          model$err.rate[,"0"],
          model$err.rate[,"1"]))

randf_error <- ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type)) + 
  theme_bw()

pdf("20240528/res/fig5/Machine_learning/randf_error.pdf", height = 8, width = 10)
print(randf_error)
dev.off()

#确定mtry值
oob.values <- vector(length=20)
for(i in 1:20) {
  temp.model <- randomForest(group ~ ., data = infr, mtry=i, ntree=1000)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate),1]
}
oob.values

# 最终确定的结果
model <- randomForest(group~.,data = infr, mtry = 3, ntree = 600,
                      proximity = TRUE, importance = TRUE)
rfi1 <- importance(model, type = 1)
rfi2 <- importance(model, type = 2)

rfi1 <- rfi1[order(rfi1, decreasing = T),]
rfi2 <- rfi2[order(rfi2, decreasing = T),]

rfif <- intersect(names(rfi1[1:top_number]), names(rfi2[1:top_number]))
pdf ("20240528/res/fig5/Machine_learning/randf.pdf", height = 8, width = 10)
varImpPlot(model)
dev.off()

fina_gene <- intersect(top.features[1:which(errors==min(errors))[1],]$FeatureName, rownames(factors))
fina_gene <- intersect(fina_gene, rfif)

write.table(fina_gene,'20240528/res/fig5/Machine_learning/hubgene.txt',row.names = F,col.names = F,quote = F)

####Venn####
#绘制韦恩图
library(venn)
venn_list <- list(top.features[1:which(errors==min(errors))[1],]$FeatureName, rownames(factors), rfif)
pdf("20240528/res/fig5/Machine_learning/venn.pdf", height = 4, width = 5)
venn(venn_list, snames = c("SVM-REF", "LASSO", "RF"), zcolor = c( "#F8766D", "#00BA38", "#619CCF"),
     ellipse = T, box = F, col = c( "#F8766D", "#00BA38", "#619CCF"), ilcs = 1, sncs = 1)
dev.off()

