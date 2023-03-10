#' Get graph theme
#' 
#' Get graph theme for plots created by 
#' \link[EWCE]{generate_bootstrap_plots_for_transcriptome}.
#' @return \code{ggplot2} graph theme.
#' @keywords internal
theme_graph <- function(){
    requireNamespace("ggplot2")
    ggplot2::theme_bw(base_size = 12, 
                                     base_family = "Helvetica") +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_line(linewidth = .25,
                                                     color = "grey"),
            axis.line = ggplot2::element_line(linewidth = .35, 
                                              color = "black"),
            text = ggplot2::element_text(size = 14),
            axis.title.x = ggplot2::element_text(vjust = -0.35),
            axis.title.y = ggplot2::element_text(vjust = 0.6),
            # legend.title = ggplot2::element_blank(), 
            strip.background = ggplot2::element_rect(fill="transparent")
        )  
}
