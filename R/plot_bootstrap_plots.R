#' Bootstrap plot
#'
#' Plot results of \link[EWCE]{generate_bootstrap_plots_for_transcriptome}.
#'
#' @return Null result.
#'
#' @keywords internal
#' 
plot_bootstrap_plots <- function(dat,
                                 tag,
                                 listFileName,
                                 cc,
                                 showGNameThresh,
                                 graph_theme,
                                 maxX,
                                 savePath) {
    requireNamespace("grDevices")
    requireNamespace("ggplot2")
    messager(cc,": Saving bootstrap plot.")
    basic_graph <- ggplot(dat, aes_string(x = "boot", y = "hit")) +
        geom_point(size = 1) +
        xlab("Mean Bootstrap Expression") +
        ylab("Expression in cell type (%)\n") +
        graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red")
    
    #### Plot without text #### 
    grDevices::pdf(sprintf(
        "%s/BootstrapPlots/qqplot_noText_%s___%s____%s.pdf",
        savePath, tag, listFileName, cc
    ), width = 3.5, height = 3.5)
    print(basic_graph + ggtitle(cc))
    grDevices::dev.off()

    # If a gene has over 25% of it's expression proportion in a cell type,
    # then list the genename
    dat$symLab <-
        ifelse(dat$hit > showGNameThresh, sprintf("  %s", dat$Gnames), "")

    basic_graph <- ggplot(dat, aes_string(x = "boot", y = "hit")) +
        geom_point(size = 2) +
        xlab("Mean Bootstrap Expression") +
        ylab("Expression in cell type (%)\n") +
        graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red")

    # Plot with gene names
    grDevices::pdf(sprintf(
        "%s/BootstrapPlots/qqplot_wtGSym_%s___%s____%s.pdf",
        savePath, tag, listFileName, cc
    ), width = 3.5, height = 3.5)
    print(basic_graph +
        geom_text(aes_string(label = "symLab"),
            hjust = 0, vjust = 0, size = 3
        ) + xlim(c(0, maxX)) +
        ggtitle(cc))
    grDevices::dev.off()

    # Plot BIG with gene names
    grDevices::pdf(sprintf(
        "%s/BootstrapPlots/qqplot_wtGSymBIG_%s___%s____%s.pdf",
        savePath, tag, listFileName, cc
    ), width = 15, height = 15)
    print(basic_graph +
        geom_text(aes_string(label = "symLab"),
            hjust = 0,
            vjust = 0, size = 3
        ) + xlim(c(0, maxX)) + ggtitle(cc))
    out <- grDevices::dev.off()
}
