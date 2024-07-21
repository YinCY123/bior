library(rWikiPathways)
library(magrittr)

hsa_pathways <- listPathwayIds(organism = "Homo sapiens")
hsa_pathways %>% str

genes <- "55107"

target_pathways <- vector()
for( pathway in hsa_pathways){
    try(chebi_ids <- getXrefList(pathway = pathway, systemCode = "L"), silent = T)
    if(length(intersect(chebi_ids, kegg_ids$ChEBI)) > 0){
        target_pathways <- append(target_pathways, pathway)
    }
    
    Sys.sleep(rnorm(1, 5, 1))
}

which(pathway == hsa_pathways)