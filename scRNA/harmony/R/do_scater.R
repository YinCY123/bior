colors_use <- c(`jurkat` = '#810F7C', `t293` = '#D09E2D',`half` = '#006D2C')

do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE,
                       do_labels = TRUE, nice_names, 
                       palette_use = colors_use,
                       pt_size = 4, point_size = .5, base_size = 12, 
                       do_points = TRUE, do_density = FALSE, h = 6, w = 8) {
    umap_use <- umap_use[, 1:2]
    colnames(umap_use) <- c('X1', 'X2')
    plt_df <- umap_use %>% data.frame() %>% 
        cbind(meta_data) %>% 
        dplyr::sample_frac(1L) 
    plt_df$given_name <- plt_df[[label_name]]
    
    if (!missing(nice_names)) {
        plt_df %<>%
            dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))
        
        plt_df[[label_name]] <- plt_df$nice_name        
    }
    
    plt <- plt_df %>% 
        ggplot2::ggplot(aes_string("X1", "X2", col = label_name, fill = label_name)) + 
        theme_test(base_size = base_size) + 
        theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1,
                                                        shape = 16, size = 4)), 
               alpha = FALSE) +
        scale_color_manual(values = palette_use) + 
        scale_fill_manual(values = palette_use) +    
        theme(plot.title = element_text(hjust = .5)) + 
        labs(x = "PC 1", y = "PC 2") 
    
    if (do_points) 
        plt <- plt + geom_point(shape = '.')
    if (do_density) 
        plt <- plt + geom_density_2d()    
    
    
    if (no_guides)
        plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
    
    if (do_labels) {
        data_labels <- plt_df %>% 
            dplyr::group_by_(label_name) %>% 
            dplyr::summarise(X1 = mean(X1), X2 = mean(X2)) %>% 
            dplyr::ungroup()
        
        plt <- plt + geom_label(data = data_labels, label.size = NA,
                                aes_string(label = label_name), 
                                color = "white", size = pt_size, alpha = 1,
                                segment.size = 0) +
            guides(col = FALSE, fill = FALSE)
    }
    
    return(plt)
}
