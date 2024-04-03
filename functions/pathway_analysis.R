library(AnnotationHub)
library(KEGGREST)
library(pathview)
library(magrittr)
library(stringr)


ah <- AnnotationHub()
# query(ah, "metabolite")
df <- ah[["AH83115"]]

# reads in compounds
compounds <- read.table("/mnt/c/Users/yincy/Downloads/compounds.csv", header = F) %>% 
    dplyr::pull(V1) %>% str_to_lower()

# get KEGG id
df$Name <- str_to_lower(df$Name)
kegg_ids <- df %>% 
    dplyr::filter(Name %in% compounds) %>% 
    dplyr::select(ChEBI, Name)

pathways <- keggLink("pathway", kegg_ids$KEGG) %>% 
    str_replace("^path:map", "")

values <- setNames(rep(1, length(unique(kegg_ids$KEGG))), kegg_ids$KEGG %>% unique)

pathview(cpd.data = values, 
         pathway.id = pathways[1:3], 
         species = "Homo sapiens", 
         plot.col.key = F, 
         limit = list(gene = 1, cpd = c(-1, 1)), 
         low = list(gene = "blue", cpd = "blue"), 
         mid = list(gene = "white", cpd = "white"), 
         high = list(gene = "red", cpd = "red"))


