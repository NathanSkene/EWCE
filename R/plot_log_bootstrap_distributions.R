#' Plot log bootstrap distributions
#'
#' Plot results of \link[EWCE]{generate_bootstrap_plots_for_transcriptome}.
#'
#' @return Null result.
#'
#' @keywords internal
plot_log_bootstrap_distributions <- function(dat,
                                             exp_mats,
                                             cc,
                                             hit_exp,
                                             tag,
                                             listFileName,
                                             graph_theme,
                                             save_dir = file.path(
                                                 tempdir(),
                                                 paste0("BootstrapPlots",
                                                        "_for_transcriptome")),
                                             height=3.5,
                                             width=3.5) {
    requireNamespace("ggplot2")
    requireNamespace("reshape2")
    
    messager(cc,": Saving log bootstrap plot with distributions.")
    #### Save path ####
    pdf_path <- file.path(
        save_dir,
        sprintf("bootDists_LOG_%s___%s____%s.pdf",
                tag, listFileName, cc
        ))
    dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    # - First get the ordered gene names
    rownames(dat) <- dat$Gnames
    datOrdered <- data.frame(GSym = rownames(dat), Pos = seq_len(dim(dat)[1]))
    # - Arrange the data frame for plotting
    melt_boot <- reshape2::melt(as.matrix(exp_mats[[cc]]))
    colnames(melt_boot) <- c("Rep", "Pos", "Exp")
    melt_boot$Exp <- melt_boot$Exp * 100
    melt_boot <- merge(melt_boot, datOrdered, by = "Pos")
    melt_boot$GSym <- factor(as.character(melt_boot$GSym),
        levels = as.character(datOrdered$GSym)
    ) 
    # - Prepare the values of the list genes to be plotted as red dots
    actVals <- data.frame(
        Pos = as.factor(seq_len(length(hit_exp))),
        vals = hit_exp * 100
    )
    actVals <- merge(actVals, datOrdered, by = "Pos")
    actVals$GSym <- factor(as.character(actVals$GSym),
        levels = as.character(datOrdered$GSym)
    ) 
    # - Determine whether changes are significant
    p <- rep(1, max(melt_boot$Pos))
    for (i in seq_len(max(melt_boot$Pos))) {
        p[i] <- sum(actVals[actVals$Pos == i, "vals"] <
            melt_boot[melt_boot$Pos == i, "Exp"]) /
            length(melt_boot[melt_boot$Pos == i, "Exp"])
    }
    ast <- rep("*", max(melt_boot$Pos))
    ast[p > 0.05] <- ""
    actVals <- cbind(actVals[order(actVals$Pos), ], ast)
    # - Plot the graph!
    wd <- 1 + length(unique(melt_boot[, 4])) * 0.175
    
    gp <-  ggplot(melt_boot) +
        geom_boxplot(aes(x = .data$GSym, y = .data$Exp), outlier.size = 0) +
        graph_theme +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_point(aes(x = .data$GSym, y = .data$vals),
                   col = "red",
                   data = actVals
        ) +
        geom_text(aes(x = .data$GSym, y = .data$vals, label = .data$ast),
                  colour = "black", col = "black", data = actVals
        ) +
        ylab("Expression in cell type (%)\n") +
        xlab("Least specific --> Most specific") +
        scale_y_log10(labels = myScalesComma, limits = c(0.01, 100))
    return(list(plot=gp,
                path=pdf_path))
}
