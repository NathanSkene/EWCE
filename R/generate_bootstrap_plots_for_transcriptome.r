#' Generate bootstrap plots
#'
#' \code{generate_bootstrap_plots_for_transcriptome} takes a genelist and a
#' single cell type transcriptome dataset and generates plots which show how
#' the expression of the genes in the list compares to those in randomly
#' generated gene lists
#'
#' @param full_results The full output of
#' \link[EWCE]{ewce_expression_data} for the same gene list.
#' @param listFileName String used as the root for files saved using
#' this function.
#' @param tt Differential expression table.
#' Can be output of \link[limma]{topTable} function.
#'  Minimum requirement is that one column stores a metric of
#' increased/decreased expression (i.e. log fold change, t-statistic for
#' differential expression etc) and another contains gene symbols.
#' @param thresh The number of up- and down- regulated genes to be included in
#' each analysis (Default: 250)
#' @param annotLevel an integer indicating which level of the annotation to
#' analyse (Default: 1).
#' @param showGNameThresh Integer. If a gene has over X percent of it's
#' expression proportion in a cell type, then list the gene name
#' @param ttSpecies Either 'mouse' or 'human' depending on which species the
#' differential expression table was generated from
#' @param sctSpecies Either 'mouse' or 'human' depending on which species the
#' single cell data was generated from.
#' @param sortBy Column name of metric in tt which should be used to sort up-
#' from down- regulated genes (Default: "t").
#' @param onlySignif Should plots only be generated for cells which have
#' significant changes?
#' @param savePath Directory where the \emph{BootstrapPlots} folder
#'  should be saved, default is a temp directory.
#' @inheritParams bootstrap_enrichment_test
#'
#' @return Saves a set of PDF files containing graphs and returns the file where
#' they are saved. These will be saved with the filename adjusted using the
#' value of \code{listFileName.} The files are saved into the
#' \emph{BootstrapPlot} folder.
#' Files start with one of the following:
#' \itemize{
#'   \item \code{qqplot_noText}: sorts the gene list according to how enriched
#'   it is in the relevant cell type. Plots the value in the target list against
#'   the mean value in the bootstrapped lists.
#'   \item \code{qqplot_wtGSym}: as above but labels the gene symbols for the
#'   highest expressed genes.
#'   \item \code{bootDists}: rather than just showing the mean of the
#'   bootstrapped lists, a boxplot shows the distribution of values
#'   \item \code{bootDists_LOG}: shows the bootstrapped distributions with the
#'   y-axis shown on a log scale
#' }
#'
#' @examples
#' \dontrun{
#' ## Load the single cell data
#' ctd <- ewceData::ctd()
#'
#' ## Set the parameters for the analysis
#' ## Use 3 bootstrap lists for speed, for publishable analysis use >10,000
#' reps <- 3
#' annotLevel <- 1 # <- Use cell level annotations (i.e. Interneurons)
#' ## Use 5 up/down regulated genes (thresh) for speed, default is 250
#' thresh <- 5
#'
#' ## Load the top table
#' tt_alzh <- ewceData::tt_alzh()
#'
#' tt_results <- EWCE::ewce_expression_data(
#'     sct_data = ctd,
#'     tt = tt_alzh,
#'     annotLevel = 1,
#'     thresh = thresh,
#'     reps = reps,
#'     ttSpecies = "human",
#'     sctSpecies = "mouse"
#' )
#'
#' ## Bootstrap significance test,
#' ## no control for transcript length or GC content
#' full_results <- EWCE::generate_bootstrap_plots_for_transcriptome(
#'     sct_data = ctd,
#'     tt = tt_alzh,
#'     thresh = thresh,
#'     annotLevel = 1,
#'     full_results = tt_results,
#'     listFileName = "examples",
#'     reps = reps,
#'     ttSpecies = "human",
#'     sctSpecies = "mouse",
#'     savePath = tempdir()
#' )
#' }
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom scales comma
generate_bootstrap_plots_for_transcriptome <- function(sct_data,
                                                       tt,
                                                       thresh = 250,
                                                       annotLevel = 1,
                                                       reps,
                                                       full_results = NA,
                                                       listFileName = "",
                                                       showGNameThresh = 25,
                                                       ttSpecies = "mouse",
                                                       sctSpecies = "mouse",
                                                       sortBy = "t",
                                                       onlySignif = TRUE,
                                                       savePath = tempdir()) {
    requireNamespace("grDevices")
    tt <- check_args_for_bootstrap_plot_generation(
        sct_data = sct_data,
        tt = tt,
        thresh = thresh,
        annotLevel = annotLevel,
        reps = reps,
        full_results = full_results,
        listFileName = listFileName,
        showGNameThresh = showGNameThresh,
        ttSpecies = ttSpecies,
        sctSpecies = sctSpecies,
        sortBy = sortBy
    )

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
                mouse.hits %in% rownames(sct_data[[annotLevel]]$specificity)
            ]
        mouse.bg <- as.character(unique(tt$MGI.symbol))
        mouse.bg <- mouse.bg[!mouse.bg %in% mouse.hits]
        mouse.bg <- mouse.bg[
            mouse.bg %in% rownames(sct_data[[annotLevel]]$specificity)
        ]
        combinedGenes <- unique(c(mouse.hits, mouse.bg))

        # Get expression data of bootstrapped genes
        if (onlySignif) {
            signif_res <- as.character(results$CellType)[results$p < 0.05]
        } else {
            signif_res <- as.character(results$CellType)
        }
        exp_mats <- get_exp_data_for_bootstrapped_genes(
            results = results,
            signif_res = signif_res,
            sct_data = sct_data,
            mouse.hits = mouse.hits,
            combinedGenes = combinedGenes,
            annotLevel = annotLevel,
            nReps = reps
        )

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

        if (!file.exists(sprintf("%s/BootstrapPlots", savePath))) {
            dir.create(file.path(savePath, "BootstrapPlots"),
                showWarnings = FALSE, recursive = TRUE
            )
        }

        tag <- sprintf("thresh%s__dir%s", thresh, dirS)

        # Plot the QQ plots
        for (cc in signif_res) {
            mean_boot_exp <- apply(exp_mats[[cc]], 2, mean)
            hit_exp <- sort(hit.exp[, cc])
            hit_exp_names <- rownames(hit.exp)[order(hit.exp[, cc])]
            dat <- data.frame(
                boot = mean_boot_exp, hit = hit_exp,
                Gnames = hit_exp_names
            )
            dat$hit <- dat$hit * 100
            dat$boot <- dat$boot * 100
            maxHit <- max(dat$hit, na.rm = TRUE)
            maxX <- max(dat$boot, na.rm = TRUE) +
                0.1 * max(dat$boot, na.rm = TRUE)

            # Plot several variants of the graph
            basic_graph <- plot_bootstrap_plots(
                dat = dat,
                tag = tag,
                listFileName = listFileName,
                cc = cc,
                showGNameThresh = showGNameThresh,
                graph_theme = graph_theme,
                maxX = maxX,
                savePath = savePath
            )

            # Plot with bootstrap distribution
            plot_with_bootstrap_distributions(
                exp_mats = exp_mats,
                cc = cc,
                hit_exp = hit_exp,
                tag = tag,
                listFileName = listFileName,
                graph_theme = graph_theme,
                savePath = savePath
            )

            # Plot with LOG bootstrap distribution
            plot_log_bootstrap_distributions(
                dat = dat,
                exp_mats = exp_mats,
                cc = cc,
                hit_exp = hit_exp,
                tag = tag,
                listFileName = listFileName,
                graph_theme = graph_theme,
                savePath = savePath
            )
        }
    }
    # return path to the saved directory in case tempdir() used
    return(savePath)
}
