####TMB####
if(!file.exists(maf_path)){
  print("maf_path不存在，自动下载")
  library(TCGAbiolinks)
  query_exp = GDCquery(project = project,
                       data.category = "Simple Nucleotide Variation",
                       data.type="Masked Somatic Mutation")
  
  GDCdownload(query = query_exp,directory = maf_path)
}

dir.create('TMB')
dir.create("TMB/plot")
dir.create("TMB/data")
dir.create("TMB/Table")
filepath = dir(path =paste0(maf_path,'/',project,'/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/'),
               pattern = "maf.gz$",
               full.names = T,
               recursive = T)

x = lapply(filepath, data.table::fread, skip = "Hugo_Symbol")
x = data.table::rbindlist(l = x, use.names = TRUE, fill = TRUE)
luad = maftools::read.maf(maf = x)

luad # 直接输入MAF对象可以查看MAF文件的基本信息
getSampleSummary(luad) # 显示样品的统计
getGeneSummary(luad) # 显示基因的统计
getClinicalData(luad) # 显示样品关联的临床数据
getFields(luad) # 显示MAF文件中的所有字段
write.mafSummary(maf=luad, basename=paste0("TMB/data/luad")) # 将详细统计结果输出到文件
while (!is.null(dev.list()))  dev.off()
pdf('TMB/plot/mafSummary.pdf',width = 10,height = 10)
plotmafSummary(maf=luad, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE)
dev.off()

maf<-read.maf(paste0('TMB/data/luad_maftools.maf'))
dat<-maf@data
sample1=unique(as.character(dat$Tumor_Sample_Barcode))
sample2=substr(sample1,1,12)
sample1=sample1[!duplicated(sample2)]
sample2=sample2[!duplicated(sample2)]
hname=names(group)[group==case]
lname=names(group)[group==control]

samp=intersect(sample2,hname)
name1=sample1[which(sample2 %in% samp)]
if(length(name1)==0)stop("case组没有SNV信息，请检查数据")
samp=intersect(sample2,lname)
name2=sample1[which(sample2 %in% samp)]
if(length(name2)==0)stop("control组没有SNV信息，请检查数据")

high_maf=subsetMaf(maf,tsb = name1)
low_maf=subsetMaf(maf,tsb = name2)
maf1=subsetMaf(maf,tsb = c(name1,name2))

save(maf1,file = "TMB/data/maf1.rda")
save(high_maf,low_maf,file="TMB/data/maf_group.rda")

a<-maf1@data
b<-unique(a[,c('Hugo_Symbol','Tumor_Sample_Barcode')])
times<-as.data.frame(table(b$Hugo_Symbol))
times<-times[order(times$Freq,decreasing = T),]
top20_gene<-as.character(times$Var1[1:20])
while (!is.null(dev.list()))  dev.off()
pdf(paste0('TMB/plot/','all','_maf.pdf'),width = 7,height = 7)
oncoplot(maf=maf1, borderCol=NULL,gene_mar = 8,genes = top20_gene,keepGeneOrder = T)
dev.off()
while (!is.null(dev.list()))  dev.off()
pdf(paste0('TMB/plot/',case,'_maf.pdf'),width = 7,height = 7)
oncoplot(maf=high_maf, borderCol=NULL,gene_mar = 8,genes = top20_gene,keepGeneOrder = T)
dev.off()
while (!is.null(dev.list()))  dev.off()
pdf(paste0('TMB/plot/',control,'_maf.pdf'),width = 7,height = 7)
oncoplot(maf=low_maf, borderCol=NULL,gene_mar = 8,genes = top20_gene,keepGeneOrder = T)
dev.off()

#肿瘤突变负荷计算
#maftools版本2.10.05
tmb1<-tmb(maf = maf1)

tmb1$risk<-ifelse(tmb1$Tumor_Sample_Barcode%in%name1,case,control)
write.csv(tmb1,'TMB/Table/tmb.csv')

tmbplot<-ggplot(data = tmb1, aes(x = risk, y = total_perMB_log, fill = risk))+ 
  #scale_fill_manual(values = c("class1"="#02786A", "class2"="#EB4B17")) +
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2,
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("log10 (Tumor mutation burden)") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  stat_compare_means(aes(group = risk,label = ..p.signif..),bracket.size = 0.6, size = 3,
                     label.x = 1.5,
                     label.y = max(tmb1$total_perMB_log)*0.95+min(tmb1$total_perMB_log)*0.05, #p值位置
                     hide.ns = T,#隐藏ns
                     method = ifelse(length(unique(tmb1$risk))==2,"wilcox.test","kruskal.test")) 
while (!is.null(dev.list()))  dev.off()
pdf("TMB/plot/TMB_group_boxplot.pdf",width = 6,height = 6)
print(tmbplot)
dev.off()

