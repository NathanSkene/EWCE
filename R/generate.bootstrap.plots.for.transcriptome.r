#' Generate bootstrap plots
#'
#' \code{generate.bootstrap.plots.for.transcriptome} takes a genelist and a 
#' single cell type transcriptome dataset and generates plots which show how 
#' the expression of the genes in the list compares to those in randomly 
#' generated gene lists
#'
#' @param sct_data List generated using \code{\link{generate.celltype.data}}
#' @param reps Number of random gene lists to generate (default=100 but should 
#' be over 10000 for publication quality results)
#' @param full_results The full output of 
#' \code{\link{bootstrap.enrichment.test}} for the same genelist
#' @param listFileName String used as the root for files saved using 
#' this function
#' @param tt Differential expression table. Can be output of limma::topTable 
#' function. Minimum requirement is that one column stores a metric of 
#' increased/decreased expression (i.e. log fold change, t-statistic for 
#' differential expression etc) and another contains either HGNC or MGI symbols.
#' @param thresh The number of up- and down- regulated genes to be included in 
#' each analysis. Dafault=250
#' @param annotLevel an integer indicating which level of the annotation to 
#' analyse. Default = 1.
#' @param showGNameThresh Integer. If a gene has over X percent of it's 
#' expression proportion in a celltype, then list the gene name
#' @param ttSpecies Either 'mouse' or 'human' depending on which species the 
#' differential expression table was generated from
#' @param sctSpecies Either 'mouse' or 'human' depending on which species the 
#' single cell data was generated from
#' @param sortBy Column name of metric in tt which should be used to sort up- 
#' from down- regulated genes. Default="t"
#' @param onlySignif Should plots only be generated for cells which have 
#' significant changes?
#' @return Saves a set of pdf files containing graphs. These will be saved 
#' with the filename adjusted using the value of listFileName. The files are 
#' saved into the 'BootstrapPlot' folder. Files start with one of the following:
#' \itemize{
#'   \item \code{qqplot_noText}: sorts the gene list according to how enriched 
#'   it is in the relevant celltype. Plots the value in the target list against 
#'   the mean value in the bootstrapped lists.
#'   \item \code{qqplot_wtGSym}: as above but labels the gene symbols for the 
#'   highest expressed genes.
#'   \item \code{bootDists}: rather than just showing the mean of the 
#'   bootstrapped lists, a boxplot shows the distribution of values
#'   \item \code{bootDists_LOG}: shows the bootstrapped distributions with the 
#'   y-axis shown on a log scale
#' }
#'
#'
#' @examples
#' # Load the single cell data
#' ctd <- ctd()
#'
#' # Set the parameters for the analysis
#' # Use 100 bootstrap lists for speed, for publishable analysis use >10000
#' reps <- 100 
#' annotLevel <- 1 # <- Use cell level annotations (i.e. Interneurons)
#'
#' # Load the top table
#' tt_alzh <- tt_alzh()
#'
#' tt_results <- ewce_expression_data(
#'     sct_data = ctd, tt = tt_alzh, annotLevel = 1,
#'     ttSpecies = "human", sctSpecies = "mouse"
#' )
#'
#' # Bootstrap significance test, no control for transcript length or GC content
#' full_results <- generate.bootstrap.plots.for.transcriptome(
#'     sct_data = ctd, tt = tt_alzh, annotLevel = 1,
#'     full_results = tt_results, listFileName = "examples", reps = reps, 
#'     ttSpecies = "human", sctSpecies = "mouse"
#' )
#' @export
#' @import ggplot2
#' @import ewceData
#' @importFrom reshape2 melt
#' @importFrom scales comma
# @import plyr
#' @import grDevices
generate.bootstrap.plots.for.transcriptome <- function(sct_data, tt, 
                                                        thresh = 250, 
                                                        annotLevel = 1, 
                                                        reps, full_results = NA,
                                                        listFileName = "", 
                                                        showGNameThresh = 25, 
                                                        ttSpecies = "mouse", 
                                                        sctSpecies = "mouse", 
                                                        sortBy = "t", 
                                                        onlySignif = TRUE) {
    tt <- check_args_for_bootstrap_plot_generation(sct_data, tt, thresh, 
                                                    annotLevel, reps, 
                                                    full_results, listFileName, 
                                                    showGNameThresh, ttSpecies, 
                                                    sctSpecies, sortBy)

    for (dirS in c("Up", "Down")) {
        a <- full_results$joint_results
        results <- a[as.character(a$Direction) == dirS, ]

        # Drop genes lacking expression data
        if (dirS == "Up") {
            tt <- tt[order(tt[, sortBy], decreasing = TRUE), ]
        }
        if (dirS == "Down") {
            tt <- tt[order(tt[, sortBy], decreasing = FALSE), ]
        }
        mouse.hits <- as.character(unique(tt$MGI.symbol[seq_len(thresh)]))
        mouse.hits <- 
            mouse.hits[
                mouse.hits %in% rownames(sct_data[[annotLevel]]$specificity)]
        mouse.bg <- as.character(unique(tt$MGI.symbol))
        mouse.bg <- mouse.bg[!mouse.bg %in% mouse.hits]
        mouse.bg <- mouse.bg[
            mouse.bg %in% rownames(sct_data[[annotLevel]]$specificity)]
        combinedGenes <- unique(c(mouse.hits, mouse.bg))

        # Get expression data of bootstrapped genes
        if (onlySignif) {
            signif_res <- as.character(results$CellType)[results$p < 0.05]
        } else {
            signif_res <- as.character(results$CellType)
        }
        exp_mats <- get_exp_data_for_bootstrapped_genes(results, signif_res, 
                                                            sct_data, 
                                                            mouse.hits, 
                                                            combinedGenes, 
                                                            annotLevel, 
                                                            nReps = reps)

        # Get expression levels of the hit genes
        hit.exp <- sct_data[[annotLevel]]$specificity[mouse.hits, ] 

        graph_theme <- theme_bw(base_size = 12, base_family = "Helvetica") +
            theme(
                panel.grid.major = element_line(size = .5, color = "grey"),
                axis.line = element_line(size = .7, color = "black"), 
                legend.position = c(0.75, 0.7), text = element_text(size = 14),
                axis.title.x = element_text(vjust = -0.35), 
                axis.title.y = element_text(vjust = 0.6)
            ) + theme(legend.title = element_blank())

        if (!file.exists("BootstrapPlots")) {
            dir.create(file.path(getwd(), "BootstrapPlots"))
        }

        tag <- sprintf("thresh%s__dir%s", thresh, dirS)

        # Plot the QQ plots
        for (cc in signif_res) {
            mean_boot_exp <- apply(exp_mats[[cc]], 2, mean)
            hit_exp <- sort(hit.exp[, cc])
            hit_exp_names <- rownames(hit.exp)[order(hit.exp[, cc])] 
            dat <- data.frame(boot = mean_boot_exp, hit = hit_exp, 
                                Gnames = hit_exp_names)
            dat$hit <- dat$hit * 100
            dat$boot <- dat$boot * 100
            maxHit <- max(dat$hit)
            maxX <- max(dat$boot) + 0.1 * max(dat$boot)

            # Plot several variants of the graph
            basic_graph <- plot_bootstrap_plots(dat, tag, listFileName, cc, 
                                                    showGNameThresh, 
                                                    graph_theme, maxX)

            # Plot with bootstrap distribution
            plot_with_bootstrap_distributions(exp_mats, cc, hit_exp, tag, 
                                                listFileName, graph_theme)

            # Plot with LOG bootstrap distribution
            plot_log_bootstrap_distributions(dat, exp_mats, cc, hit_exp, tag, 
                                                listFileName, graph_theme)
        }
    }
}

check_args_for_bootstrap_plot_generation <- function(sct_data, tt, thresh, 
                                                        annotLevel, reps, 
                                                        full_results, 
                                                        listFileName, 
                                                        showGNameThresh, 
                                                        ttSpecies, sctSpecies, 
                                                        sortBy) {
    # Check the arguments
    correct_length <- length(full_results) == 5
    required_names <- c("joint_results", "hit.cells.up", 
                            "hit.cells.down", "bootstrap_data.up", 
                            "bootstrap_data.down")
    all_required_names <- sum(names(full_results) %in% required_names) == 5
    err_msg <- paste0("ERROR: full_results is not valid output from the",
                        " ewce_expression_data function. This function only",
                        " takes data generated from transcriptome analyses.")
    if (!correct_length | !all_required_names) {
        stop(err_msg)
    }

    # Check the arguments
    err_msg2 <- paste0("ERROR: tt does not contain a column with value",
                        " passed in sortBy argument")
    if (!sortBy %in% colnames(tt)) {
        stop(err_msg2)
    }
    err_msg3 <- paste0("ERROR: length of table is less than twice the",
                        " size of threshold")
    if (dim(tt)[1] < (thresh * 2)) {
        stop(err_msg3)
    }

    # Check that the top table has correct columns
    err_msg4 <- paste0("ERROR: if ttSpecies==human then there must be an",
                        " HGNC.symbol column")
    err_msg5 <- paste0("ERROR: if ttSpecies==human then there must be an",
                        " MGI.symbol column")
    if (ttSpecies == "human" & !"HGNC.symbol" %in% colnames(tt)) {
        stop(err_msg4)
    }
    if (ttSpecies == "mouse" & !"MGI.symbol" %in% colnames(tt)) {
        stop(err_msg5)
    }

    if (ttSpecies == "human" & sctSpecies == "human") {
        tt$MGI.symbol <- tt$HGNC.symbol
    }
    m2h <- ewceData::mouse_to_human_homologs()[, c("MGI.symbol", "HGNC.symbol")]
    if (ttSpecies == "human" & sctSpecies == "mouse") {
        tt <- merge(tt, m2h, by = "HGNC.symbol")
    }
    if (ttSpecies == "mouse" & sctSpecies == "human") {
        tt <- merge(tt, m2h, by = "HGNC.symbol")
        tt$MGI.symbol <- tt$HGNC.symbol
    }
    return(tt)
}


get_exp_data_for_bootstrapped_genes <- function(results, signif_res, sct_data, 
                                                    mouse.hits, combinedGenes, 
                                                    annotLevel, nReps = NA) {
    exp_mats <- list()
    for (cc in signif_res) {
        exp_mats[[cc]] <- matrix(0, nrow = nReps, ncol = length(mouse.hits))
        rownames(exp_mats[[cc]]) <- sprintf("Rep%s", seq_len(nReps))
    }
    for (s in seq_len(nReps)) {
        bootstrap_set <- sample(combinedGenes, length(mouse.hits))
        ValidGenes <- rownames(sct_data[[annotLevel]]$specificity)[
            rownames(sct_data[[annotLevel]]$specificity) %in% bootstrap_set]

        expD <- sct_data[[annotLevel]]$specificity[ValidGenes, ]

        for (cc in signif_res) {
            exp_mats[[cc]][s, ] <- sort(expD[, cc])
        }
    }
    return(exp_mats)
}

plot_with_bootstrap_distributions <- function(exp_mats, cc, hit_exp, tag, 
                                                listFileName, graph_theme) {
    melt_boot <- reshape2::melt(exp_mats[[cc]])
    colnames(melt_boot) <- c("Rep", "Pos", "Exp")
    actVals <- data.frame(pos = as.factor(seq_len(length(hit_exp))), 
                            vals = hit_exp)
    pdf(sprintf("BootstrapPlots/bootDists_%s___%s____%s.pdf", 
                    tag, listFileName, cc), width = 3.5, height = 3.5)
    melt_boot$Pos <- as.factor(melt_boot$Pos)
    print(ggplot(melt_boot) +
        geom_boxplot(aes_string(x = "Pos", y = "Exp"), outlier.size = 0) +
        geom_point(aes_string(x = "pos", y = "vals"), 
                    col = "red", data = actVals) +
        ylab("Expression in cell type (%)\n") +
        xlab("Least specific --> Most specific") +
        scale_x_discrete(breaks = NULL) +
        graph_theme)
    dev.off()
}

plot_bootstrap_plots <- function(dat, tag, listFileName, cc, showGNameThresh, 
                                    graph_theme, maxX) {
    basic_graph <- ggplot(dat, aes_string(x = "boot", y = "hit")) +
        geom_point(size = 1) +
        xlab("Mean Bootstrap Expression") +
        ylab("Expression in cell type (%)\n") +
        graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red")

    # Plot without text
    pdf(sprintf("BootstrapPlots/qqplot_noText_%s____%s____%s.pdf", 
                    tag, listFileName, cc), width = 3.5, height = 3.5)
    print(basic_graph + ggtitle(cc))
    dev.off()

    # If a gene has over 25% of it's expression proportion in a celltype, 
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
    pdf(sprintf("BootstrapPlots/qqplot_wtGSym_%s___%s____%s.pdf", 
                    tag, listFileName, cc), width = 3.5, height = 3.5)
    print(basic_graph +
        geom_text(aes_string(label = "symLab"), 
                    hjust = 0, vjust = 0, size = 3) + xlim(c(0, maxX)) + 
            ggtitle(cc))
    dev.off()

    # Plot BIG with gene names
    pdf(sprintf("BootstrapPlots/qqplot_wtGSymBIG_%s___%s____%s.pdf", tag, 
                    listFileName, cc), width = 15, height = 15)
    print(basic_graph +
        geom_text(aes_string(label = "symLab"), hjust = 0, 
                    vjust = 0, size = 3) + xlim(c(0, maxX)) + ggtitle(cc))
    dev.off()

    return(basic_graph)
}

myScalesComma <- function(x) {
    return(scales::comma(x = x, accuracy = 0.01))
}

plot_log_bootstrap_distributions <- function(dat, exp_mats, cc, hit_exp, tag, 
                                                listFileName, graph_theme) {
    # - First get the ordered gene names
    rownames(dat) <- dat$Gnames
    datOrdered <- data.frame(GSym = rownames(dat), Pos = seq_len(dim(dat)[1]))

    # - Arrange the data frame for plotting
    melt_boot <- reshape2::melt(exp_mats[[cc]])
    colnames(melt_boot) <- c("Rep", "Pos", "Exp")
    melt_boot$Exp <- melt_boot$Exp * 100
    melt_boot <- merge(melt_boot, datOrdered, by = "Pos")
    melt_boot$GSym <- factor(as.character(melt_boot$GSym), 
                                levels = as.character(datOrdered$GSym))

    # - Prepare the values of the list genes to be plotted as red dots
    actVals <- data.frame(Pos = as.factor(seq_len(length(hit_exp))), 
                            vals = hit_exp * 100)
    actVals <- merge(actVals, datOrdered, by = "Pos")
    actVals$GSym <- factor(as.character(actVals$GSym), 
                            levels = as.character(datOrdered$GSym))

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
    pdf(sprintf("BootstrapPlots/bootDists_LOG_%s___%s____%s.pdf", 
                    tag, listFileName, cc), width = wd, height = 4)
    print(
        ggplot(melt_boot) +
            geom_boxplot(aes_string(x = "GSym", y = "Exp"), outlier.size = 0) +
            graph_theme +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            geom_point(aes_string(x = "GSym", y = "vals"), col = "red", 
                        data = actVals) +
            geom_text(aes_string(x = "GSym", y = "vals", label = "ast"), 
                        colour = "black", col = "black", data = actVals) +
            ylab("Expression in cell type (%)\n") +
            xlab("Least specific --> Most specific") +
            scale_y_log10(labels = myScalesComma, limits = c(0.01, 100))
    )
    dev.off()
}
