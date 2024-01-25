palign <- function(pattern_dir, sub_dir, type = "global", align_out_dir){
    require(Biostrings)
    require(stringr)
    
    pattern <- read.table(pattern_dir) %>% .[, 1, drop = TRUE] %>% str_c(collapse = "") %>% DNAString() %>% DNAStringSet()
    sub <- read.table(sub_dir) %>% .[, 1, drop = TRUE] %>% str_c(collapse = "") %>% DNAString() %>% DNAStringSet()
    
    align <- pairwiseAlignment(pattern = pattern, subject = sub, type = type)
    
    writePairwiseAlignments(x = align, file = align_out_dir)
}