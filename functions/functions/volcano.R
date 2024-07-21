library(ggplot2)
library(ggrepel)
library(magrittr)

volcano <- function(data, x, y, p = 0.05, top = 10, threshold = 1, ann = TRUE, label = "symbol", ...){
    df <- data
    
    x <- colnames(data)[which(colnames(data) == x)]
    y <- colnames(data)[which(colnames(data) == y)]
    
    # df_to_plot <- transform(df, 
    #                         log10pval = -log10(y), 
    #                         colors = ifelse(y <= p & x >= threshold, "tomato", 
    #                                         ifelse(y <= p & x <= -threshold, "steelblue", "grey")), 
    #                         size = ifelse(abs(x) >= threshold & y <= p, 3, 0.5))
    
    f_to_plot <- df %>% 
        dplyr::mutate(log10pval = -log10(y), 
                      colors = ifelse(y <= p & x >= threshold, "tomato", 
                                      ifelse(y <= p & x <= -threshold, "steelblue", "grey")))
    
    
    data_to_ann <- rbind(df_to_plot %>% dplyr::filter(y <= p & x >= threshold) %>% dplyr::arrange(y) %>% head(top), 
                         df_to_plot %>% dplyr::filter(y <= p & x <= -threshold) %>% dplyr::arrange(y) %>% head(top))
    
    p <- df_to_plot %>% 
        ggplot(aes(x, log10pval)) +
        geom_point(aes(color = colors, size = size)) +
        geom_vline(xintercept = c(-x , x), linetype = 2, linewidth = 0.5, color = "grey") +
        geom_hline(yintercept = log10(p), linetype = 2, linewidth = 0.5, color = "grey") +
        scale_size_identity() +
        scale_x_continuous(name = "log2(FC)") +
        scale_y_continuous(name = "-log10(P Value)") +
        theme(legend.position = "none")
    
    if(ann){
        p <- p + geom_label_repel(data = data_to_ann, aes(x, log10pval, label = label))
    }
    return(p)
}


