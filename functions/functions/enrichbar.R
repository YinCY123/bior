library(fgsea)
library(ggplot2)
library(dplyr)

enrichbar <- function(res, top = 10, dir = NULL, text_size = 3, expand_L = 1, expand_R = 1, step = 2, width = 7, height = 7, ...){
    up <- res %>% as.data.frame %>% dplyr::filter(NES > 0) %>% dplyr::arrange(pval) %>% head(top)
    up$group <- "up"
    down <- res %>% as.data.frame %>% dplyr::filter(NES < 0) %>% dplyr::arrange(pval) %>% head(top)
    down$group <- "down"
    
    df_to_plot <- rbind(up, down) %>% 
        dplyr::mutate(log10pval = ifelse(group == "up", -log10(pval), log10(pval)), 
                      colors = ifelse(group == "up", "tomato", "steelblue")) %>% 
        dplyr::arrange(desc(log10pval)) %>% 
        dplyr::mutate(y = nrow(.):1)
    
    x_min <- min(df_to_plot$log10pval)
    x_max <- max(df_to_plot$log10pval)
    
    if(is.numeric(expand_L)){
        x_min <- min(x_min, expand_L)
    }
    
    if(is.numeric(expand_R)){
        x_max <- max(x_max, expand_R)
    }
    
    p <- df_to_plot %>% 
        ggplot(aes(log10pval, reorder(pathway, log10pval))) +
        geom_bar(aes(fill = colors), stat = "identity") +
        geom_text(data = df_to_plot[df_to_plot$group == "up", ], aes(-0.1, y, label = pathway, hjust = 1), size = text_size) +
        geom_text(data = df_to_plot[df_to_plot$group == "down", ], aes(0.1, y, label = pathway, hjust = 0), size = text_size) +
        scale_fill_identity() +
        scale_y_discrete(name = NULL, label = NULL) +
        scale_x_continuous(name = "-log10 P Value", 
                           breaks = c(rev(seq(0, x_min, -step)), seq(step, x_max, step)), 
                           labels = c(abs(rev(seq(0, x_min, -step))), seq(step, x_max, step)), 
                           limits = c(x_min, x_max)) +
        theme(axis.ticks.y = element_blank())
    if(!is.null(dir)){
        ggsave(dir)
    }else{
        return(p)
    }
}