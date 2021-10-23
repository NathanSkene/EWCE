plot_with_bootstrap_distributions <- function(exp_mats,
                                              cc,
                                              hit_exp,
                                              tag,
                                              listFileName,
                                              graph_theme,
                                              savePath) {
    requireNamespace("grDevices")
    melt_boot <- reshape2::melt(exp_mats[[cc]])
    colnames(melt_boot) <- c("Rep", "Pos", "Exp")
    actVals <- data.frame(
        pos = as.factor(seq_len(length(hit_exp))),
        vals = hit_exp
    )
    grDevices::pdf(sprintf(
        "%s/BootstrapPlots/bootDists_%s___%s____%s.pdf",
        savePath, tag, listFileName, cc
    ), width = 3.5, height = 3.5)
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
    grDevices::dev.off()
}
