#' Bootstrap plot
#'
#' Plot results of \link[EWCE]{generate_bootstrap_plots_for_transcriptome}.
#'
#' @return Null result.
#'
#' @keywords internal
#' 
bootstrap_plots_for_transcriptome <- function(dat,
                                              tag,
                                              listFileName,
                                              cc,
                                              showGNameThresh,
                                              graph_theme,
                                              maxX,
                                              save_dir = file.path(
                                                  tempdir(),
                                                  paste0("BootstrapPlots",
                                                         "_for_transcriptome")),
                                              height = 3.5,
                                              width = 3.5,
                                              show_plot = TRUE) { 
    requireNamespace("ggplot2") 
    
    # If a gene has over 25% of it's expression proportion in a cell type,
    # then show the gene name.
    dat$symLab <-ifelse(dat$hit > showGNameThresh,
                        sprintf("  %s", dat$Gnames),
                        "")
    plots <- list()
    #### Setup file paths ####
    files <- file.path(save_dir,
                       sprintf(c("qqplot_noText_%s___%s____%s.pdf",
                                 "qqplot_wtGSym_%s___%s____%s.pdf",
                                 "qqplot_wtGSymBIG_%s___%s____%s.pdf"),
                               tag, listFileName, cc)
    ) 
    for(f in files){
        dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)   
    } 
    #### Plot 1: without text ####
    basic_graph <- ggplot(dat, aes(x = .data$boot, y = .data$hit)) +
        geom_point(size = 1) +
        labs(title = cc,
             x="Mean Bootstrap Expression", 
             y="Expression in cell type (%)\n") +
        graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red")
    plots[["plot1"]] <- basic_graph
    messager(cc,": Saving bootstrap plot -->",files[[1]]) 
    ggplot2::ggsave(filename = files[[1]],
                    plot = basic_graph,
                    width = width,
                    height = height) 

    #### Plot 2: with gene name ####
    basic_graph2 <- ggplot(dat, aes(x = .data$boot, y = .data$hit)) +
        geom_point(size = 2) +
        labs(title = cc,
             x="Mean Bootstrap Expression", 
             y="Expression in cell type (%)\n") + 
        graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red") +
        geom_text(aes(label = .data$symLab),
                  hjust = 0, vjust = 0, size = 3
        ) + xlim(c(0, maxX)) 
    plots[["plot2"]] <- basic_graph2
    messager(cc,": Saving bootstrap plot -->",files[[2]]) 
    ggplot2::ggsave(filename = files[[2]],
                    plot = basic_graph2,
                    width = width,
                    height = height)   
    #### Plot 3: BIG with gene names ####
    messager(cc,": Saving bootstrap plot -->",files[[3]]) 
    ggplot2::ggsave(filename = files[[3]],
                    plot = basic_graph2,
                    width = width*4.2,
                    height = height*4.2)   
    #### Show plot ####
    if(isTRUE(show_plot)) methods::show(plots)
    return(list(plots=plots,
                files=files))
}
