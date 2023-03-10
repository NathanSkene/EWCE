#' Plot with bootstrap distributions
#'
#' Plot results of \link[EWCE]{generate_bootstrap_plots_for_transcriptome}.
#'
#' @return Null result.
#'
#' @keywords internal
plot_with_bootstrap_distributions <- function(exp_mats,
                                              cc,
                                              hit_exp,
                                              tag,
                                              listFileName,
                                              graph_theme,
                                              save_dir) {
    requireNamespace("grDevices")
    requireNamespace("ggplot2")
    requireNamespace("reshape2")
    messager(cc,": Saving bootstrap plot with distributions.")
    melt_boot <- reshape2::melt(as.matrix(exp_mats[[cc]]))
    colnames(melt_boot) <- c("Rep", "Pos", "Exp")
    actVals <- data.frame(
        pos = as.factor(seq_len(length(hit_exp))),
        vals = hit_exp
    )
  
    pdf_path <- file.path(
        save_dir,
        sprintf("bootDists_%s___%s____%s.pdf",
                tag, listFileName, cc
    ))
    grDevices::pdf(pdf_path, width = 3.5, height = 3.5)
    
    melt_boot$Pos <- as.factor(melt_boot$Pos)
    print(ggplot(melt_boot) +
        geom_boxplot(aes_string(x = "Pos", y = "Exp"), outlier.size = 0) +
        geom_point(aes_string(x = "pos", y = "vals"),
            col = "red", data = actVals
        ) +
        ylab("Expression in cell type (%)\n") +
        xlab("Least specific --> Most specific") +
        scale_x_discrete(breaks = NULL) +
        graph_theme)
    out <- grDevices::dev.off()
}
