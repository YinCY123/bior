# 单基因GSEA -----------------------------------------------------------------

number_gsea <- 5
#GSEA 需要样本中全部的基因及其logFC（差异和非差异都需要）
library(GSEABase)
library(enrichplot)
library(cowplot)
library(msigdbr)
library(ggridges)
library(clusterProfiler)
library(limma)
library(dplyr)
library(ggplot2)
msgdC2 = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")

colors <- c("#4DBBD5", "#E64B35","#00A087", "#3C5488", "#F39B7F")
pathway_gene=cbind(data.frame(msgdC2$gs_name),data.frame(msgdC2$gene_symbol))  
colnames(pathway_gene)=c('term','gene')

#关键基因是从机器学习筛选的13个预后基因中经免疫筛选得到的3个基因作为乳腺癌的诊断标志物。见scatterplot_use

#关键基因是从机器学习筛选的47-16(基因level与预后有关)个预后基因
hubgene <- read.table("20240528/res/fig5/inter_genes.txt", header = F)[[1]]
# hubgene <- genes
exp_combat <- readRDS("20240528/data/mtx_combat.rds")
# exp_combat <- exp4
dir.create('20240528/res/fig8/hubgene')
#根据单个基因在每个样本中的表达量分配高低表达标签，大于等于中位值是高表达，否则是低表达
for(i in hubgene){
  dir.create(paste0('20240528/res/fig8/hubgene/',i))
  group <- ifelse(as.numeric(exp_combat[i, ]) >= median(as.numeric(exp_combat[i,])),'High','Low')
  group <- factor(group,levels = c('High','Low'))
  design <- model.matrix(~ 0 + group)
  rownames(design) <- colnames(exp_combat)
  colnames(design) <- levels(group)
  ## 线性建模
  fit <- lmFit(exp_combat,design)
  cont.matrix <- makeContrasts(contrasts = paste0("High",'-',"Low"), levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  ## 经验贝叶斯调整
  fit2 <- eBayes(fit2)
  ## 筛选差异基因
  dif <- topTable(fit2, coef = 1, n = Inf)
  dif <- na.omit(dif)
  dif$gene_symbol <- row.names(dif)
  ge=dif$logFC
  names(ge)=dif$gene_symbol
  ge = sort(ge, decreasing = T)
  
  gsea<-GSEA(ge, 
             TERM2GENE = pathway_gene,
             # nPerm = 1000,
             verbose = FALSE,
             by = "fgsea",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05, 
             eps = 0)
  
  res <- gsea@result
  if(nrow(res)==0){print(paste0(i,'未富集到通路'));next}
  write.csv(res,paste0('20240528/res/fig8/hubgene/',i,'/',i,'_GSEA.csv'))
  
  gsea1<-gsea%>%top_n(n=number_gsea,wt=NES)

  
  p_ridgeplot <- ridgeplot(gsea1,
                           fill = "p.adjust",
                           orderBy = "NES",
                           showCategory = number_gsea,
                           label_format = 10)
  
  pr = p_ridgeplot$data %>% 
    ggplot(aes(x = value, y = category, fill = category)) +
    geom_density_ridges(scale = 2, size = 0.3, rel_min_height = 0.00, alpha = 0.5) +
    scale_fill_cyclical(values = colors) +
    labs(caption = '', x='', y='', title = i) +
    theme(axis.text = element_text(size = 6),
          axis.ticks = element_line(linewidth = 0.5), axis.line = element_line(size = 0.5),
          panel.background = element_blank(),title = element_text(size = 7,face = 'bold'),
          panel.grid = element_blank())
  
  pdf(paste0('20240528/res/fig8/hubgene/',i,'/',i,'_GSEA_ridges.pdf'),width = 12/2.54,height = 7/2.54)
  print(pr)
  while (!is.null(dev.list()))  dev.off()
}
