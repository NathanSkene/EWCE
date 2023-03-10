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
                                 save_dir) {
    requireNamespace("grDevices")
    requireNamespace("ggplot2")
   
    basic_graph <- ggplot(dat, aes_string(x = "boot", y = "hit")) +
        geom_point(size = 1) +
        xlab("Mean Bootstrap Expression") +
        ylab("Expression in cell type (%)\n") +
        graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red")
    
    
    #### Make dir ####
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    #### Plot without text #### 
    pdf_path <- file.path(
        save_dir,
        sprintf("qqplot_noText_%s___%s____%s.pdf",
                tag, listFileName, cc
    ))
    messager(cc,": Saving bootstrap plot -->",pdf_path) 
    grDevices::pdf(pdf_path, width = 3.5, height = 3.5)
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
    pdf_path <- file.path(
        save_dir,
        sprintf("qqplot_wtGSym_%s___%s____%s.pdf",
                tag, listFileName, cc
    ))
    messager(cc,": Saving bootstrap plot -->",pdf_path) 
    grDevices::pdf(pdf_path, width = 3.5, height = 3.5)
    print(basic_graph +
        geom_text(aes_string(label = "symLab"),
            hjust = 0, vjust = 0, size = 3
        ) + xlim(c(0, maxX)) +
        ggtitle(cc))
    grDevices::dev.off()

    # Plot BIG with gene names
    pdf_path <- file.path(
        save_dir,
        sprintf("qqplot_wtGSymBIG_%s___%s____%s.pdf",
                tag, listFileName, cc
    ))
    messager(cc,": Saving bootstrap plot -->",pdf_path) 
    grDevices::pdf(pdf_path, width = 15, height = 15)
    print(basic_graph +
        geom_text(aes_string(label = "symLab"),
            hjust = 0,
            vjust = 0, size = 3
        ) + xlim(c(0, maxX)) + ggtitle(cc))
    out <- grDevices::dev.off()
}
