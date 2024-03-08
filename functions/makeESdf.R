makeESdf <- function(es, features, join_by, ...){
    mtx <- exprs(es)[features, ] %>% as.data.frame() %>% tibble::rownames_to_column("ID")
    fdata <- fData(es)
    df <- dplyr::left_join(mtx, fdata, by = join_by)
    
    return(df)
}