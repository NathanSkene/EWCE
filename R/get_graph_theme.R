get_graph_theme <- function(){
    graph_theme <- ggplot2::theme_bw(base_size = 12, 
                                     base_family = "Helvetica") +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_line(size = .5, color = "grey"),
            axis.line = ggplot2::element_line(size = .7, color = "black"),
            legend.position = c(0.75, 0.7), 
            text = ggplot2::element_text(size = 14),
            axis.title.x = ggplot2::element_text(vjust = -0.35),
            axis.title.y = ggplot2::element_text(vjust = 0.6),
            legend.title = ggplot2::element_blank()
        ) 
    return(graph_theme)
}