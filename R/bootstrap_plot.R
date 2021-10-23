bootstrap_plot <- function(exp_mats,
                           hit.exp,
                           cc,
                           savePath,
                           listFileName,
                           show_plot = TRUE) {
    requireNamespace("grDevices")
    ### Setup theme ####
    graph_theme <- theme_bw(base_size = 12, base_family = "Helvetica") +
        theme(
            panel.grid.major = element_line(size = .5, color = "grey"),
            axis.line = element_line(size = .7, color = "black"),
            legend.position = c(0.75, 0.7), text = element_text(size = 14),
            axis.title.x = element_text(vjust = -0.35),
            axis.title.y = element_text(vjust = 0.6)
        ) + theme(legend.title = element_blank())

    #### Create plotting data ####
    mean_boot_exp <- apply(exp_mats[[cc]], 2, mean)
    hit_exp <- sort(hit.exp[, cc])
    hit_exp_names <- rownames(hit.exp)[order(hit.exp[, cc])]
    dat <- data.frame(
        boot = mean_boot_exp,
        hit = hit_exp,
        Gnames = hit_exp_names
    )
    dat$hit <- dat$hit * 100
    dat$boot <- dat$boot * 100
    maxHit <- max(dat$hit, na.rm = TRUE)
    maxX <- max(dat$boot, na.rm = TRUE) + 0.1 * max(dat$boot, na.rm = TRUE)
    # Plot several variants of the graph
    basic_graph <- ggplot(dat, aes_string(x = "boot", y = "hit")) +
        geom_point(size = 1) +
        xlab("Mean Bootstrap Expression") +
        ylab("Expression in cell type (%)\n") +
        graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red")

    # Plot without text
    grDevices::pdf(sprintf(
        "%s/BootstrapPlots/qqplot_noText____%s____%s.pdf", savePath,
        listFileName, cc
    ), width = 3.5, height = 3.5)
    print(basic_graph)
    grDevices::dev.off()

    dat$symLab <- ifelse(dat$hit > 25, sprintf("  %s", dat$Gnames), "")

    basic_graph <- ggplot(dat, aes_string(x = "boot", y = "hit")) +
        geom_point(size = 2) +
        xlab("Mean Bootstrap Expression") +
        ylab("Expression in cell type (%)\n") +
        graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red")

    # Plot with gene names
    grDevices::pdf(sprintf(
        "%s/BootstrapPlots/qqplot_wtGSym____%s____%s.pdf", savePath,
        listFileName, cc
    ), width = 3.5, height = 3.5)

    print(basic_graph +
        geom_text(aes_string(label = "symLab"),
            hjust = 0, vjust = 0, size = 3
        ) + xlim(c(0, maxX)))
    grDevices::dev.off()

    # Plot with bootstrap distribution
    melt_boot <- reshape2::melt(exp_mats[[cc]])
    colnames(melt_boot) <- c("Rep", "Pos", "Exp")
    actVals <- data.frame(
        pos = as.factor(seq_len(length(hit_exp))),
        vals = hit_exp
    )

    grDevices::pdf(sprintf(
        "%s/BootstrapPlots/bootDists____%s____%s.pdf", savePath,
        listFileName, cc
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

    #### Plot with LOG bootstrap distribution ####
    # - First get the ordered gene names
    rownames(dat) <- dat$Gnames
    datOrdered <- data.frame(
        GSym = rownames(dat),
        Pos = seq_len(dim(dat)[1])
    )

    # - Arrange the data frame for plotting
    melt_boot <- reshape2::melt(exp_mats[[cc]])
    colnames(melt_boot) <- c("Rep", "Pos", "Exp")
    melt_boot$Exp <- melt_boot$Exp * 100
    melt_boot <- merge(melt_boot, datOrdered, by = "Pos")
    melt_boot$GSym <- factor(as.character(melt_boot$GSym),
        levels = as.character(datOrdered$GSym)
    )

    #### - Prepare the values of the list genes to be plotted as red dots ####
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
    for (i in seq_len(max(melt_boot$Pos, na.rm = TRUE))) {
        p[i] <- sum(actVals[actVals$Pos == i, "vals"] <
            melt_boot[melt_boot$Pos == i, "Exp"]) /
            length(melt_boot[melt_boot$Pos == i, "Exp"])
    }
    ast <- rep("*", max(melt_boot$Pos, na.rm = TRUE))
    ast[p > 0.05] <- ""
    actVals <- cbind(actVals[order(actVals$Pos), ], ast)
    # - Plot the graph!
    wd <- 1 + length(unique(melt_boot[, 4])) * 0.175
    grDevices::pdf(sprintf(
        "%s/BootstrapPlots/bootDists_LOG____%s____%s.pdf",
        savePath, listFileName, cc
    ), width = wd, height = 4)
    # melt_boot$Exp=melt_boot$Exp+0.00000001
    melt_boot <- melt_boot[melt_boot$Exp != 0, ]
    gg <- ggplot(melt_boot) +
        geom_boxplot(aes_string(x = "GSym", y = "Exp"), outlier.size = 0) +
        graph_theme +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_point(aes_string(x = "GSym", y = "vals", fill = "vals"),
            data = actVals
        ) +
        scale_fill_gradient(low = "blue", high = "red") +
        geom_text(aes_string(x = "GSym", y = "log1p(vals)", label = "ast"),
            colour = "black", data = actVals
        ) +
        ylab("Expression in cell type (%)\n") +
        xlab("Least specific --> Most specific")
    print(gg)
    grDevices::dev.off()
}
