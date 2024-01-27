plotSums <- function(sces, features, feature_type = "ensembl", species = "human", color_by = "cell_type", q, dimred = "TSNE"){
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(tidyr))
    suppressPackageStartupMessages(require(scuttle))
    suppressPackageStartupMessages(require(scater))
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(magrittr))
    
    # convert feature to ensembl gene id, only support human and mouse
    if(feature_type != "ensembl"){
        if(species == "human"){
            suppressPackageStartupMessages(require(EnsDb.Hsapiens.v86))
            features <- mapIds(EnsDb.Hsapiens.v86, keys = features, keytype = feature_type, column = "GENEID")
        }else{
            suppressPackageStartupMessages(require(EnsDb.Mmusculus.v79))
            features <- mapIds(EnsDb.Mmusculus.v79, keys = features, keytype = feature_type, column = "GENEID")
        }
    }
    
    df <- makePerCellDF(x = sces, 
                        features = features, 
                        use.dimred = dimred) %>% 
        pivot_longer(cols = starts_with("ENSG"), names_to = "ensembl", values_to = "logcounts")
    
    # stat
    df <- df %>% group_by(Barcode, cell_type) %>% mutate(sum = sum(logcounts))
    
    # cell location
    if(dimred == "TSNE"){
        cell_locations <- df %>% group_by(cell_type) %>% summarise(x = median(TSNE.1), y = median(TSNE.2))
    }else{
        cell_locations <- df %>% group_by(cell_type) %>% summarise(x = median(UMAP.1), y = median(UMAP.2))
    }
    
    # visualization
    if(dimred == "TSNE"){
        df %>% ggplot(aes(TSNE.1, TSNE.2)) +
            geom_point(aes(color = sum), size = 0.5) +
            geom_density_2d(color = "grey", linewidth = 0.2) +
            geom_text(data = cell_locations, aes(x, y, label = cell_type)) +
            scale_color_gradient2(name = "energy\nstress", 
                                  low = "blue", mid = "white", high = "red", midpoint = median(df$sum))
    }else{
        df %>% ggplot(aes(UMAP.1, UMAP.2)) +
            geom_point(aes(color = ifelse(stat == "sum", sum, quantile)), size = 0.5) +
            geom_density_2d(color = "grey", linewidth = 0.2) +
            geom_text(data = cell_locations, aes(x, y, label = cell_type)) +
            scale_color_gradient2(name = "energy\nstress", 
                                  low = "blue", mid = "white", high = "red", midpoint = median(df$sum))
        }
}