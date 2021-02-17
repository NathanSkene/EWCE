#' Plot EWCE results
#'
#' \code{ewce.plot} generates plots of EWCE enrichment results
#'
#' @param total_res results dataframe generated using 
#' \code{\link{bootstrap.enrichment.test}} or 
#' \code{\link{ewce_expression_data}} functions. Multiple results tables can be 
#' merged into one results table, as long as the 'list' column is set to 
#' distinguish them.
#' @param mtc_method method to be used for multiple testing correction. 
#' Argument is passed to \code{\link{p.adjust}}. Valid options are "holm", 
#' "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr" or "none". Default method is bonferroni.
#' @param ctd Should be provided so that the dendrogram can be taken from it 
#' and added to plots
#' @return A ggplot containing the plot
#' @examples
#' # Load the single cell data
#' data(ctd, package="ewceData")
#'
#' # Set the parameters for the analysis
#' # Use 100 bootstrap lists for speed, for publishable analysis use >10000
#' reps <- 100
#'
#' # Load the gene list and get human orthologs
#' data("example_genelist", package="ewceData")
#' data("mouse_to_human_homologs", package="ewceData")
#' m2h <- unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])
#' mouse.hits <- 
#'      unique(m2h[m2h$HGNC.symbol %in% example_genelist, "MGI.symbol"])
#' mouse.bg <- unique(m2h$MGI.symbol)
#'
#' # Bootstrap significance test, no control for transcript length or GC content
#' full_results <- bootstrap.enrichment.test(
#'     sct_data = ctd, hits = mouse.hits,
#'     bg = mouse.bg, reps = reps, annotLevel = 2, 
#'     sctSpecies = "mouse", genelistSpecies = "mouse"
#' )
#'
#' # Generate the plot
#' print(ewce.plot(full_results$results, mtc_method = "BH"))
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @import stats
#' @import ggdendro
#' @import gridExtra
#' @importFrom grid unit
# @import plyr
ewce.plot <- function(total_res, mtc_method = "bonferroni", ctd = NULL) {
    err_msg <- paste0("ERROR: Invalid mtc_method argument. Please see",
                        " '?p.adjust' for valid methods.")
    if (!mtc_method %in% c("holm", "hochberg", "hommel", 
                            "bonferroni", "BH", "BY", "fdr", "none")) {
        stop(err_msg)
    }
    multiList <- TRUE
    if (is.null(total_res$list)) {
        multiList <- FALSE
    }

    # Check if ctd is provided (if so, dendrogram is to be added)
    make_dendro <- FALSE
    if (!is.null(ctd)) {
        make_dendro <- TRUE
        # If using dendrogram... Find the relevant level of the CTD annotation
        cells.in.ctd <- function(ctdIN, cells) {
            if (sum(!cells %in% colnames(ctdIN$specificity) == 0)) {
                return(1)
            } else {
                return(0)
            }
        }
        if (length(ctd[[1]]$plotting) > 0) {
            annotLevel <- 
                which(unlist(lapply(ctd, FUN = cells.in.ctd, 
                                        cells = as.character(
                                                    total_res$CellType))) == 1)
            err_msg2 <- paste0("All of the cells within total_res should come",
                                " from a single annotation layer of the CTD")
            if (length(annotLevel) == 0) {
                stop(err_msg2)
            }
        }

        # Set order of cells
        if (length(ctd[[annotLevel]]$plotting) > 0) {
            total_res$CellType <- 
                factor(total_res$CellType, 
                        levels = ctd[[annotLevel]]$plotting$cell_ordering)
        }
    }

    # Multiple testing correction across all rows
    total_res$q <- p.adjust(total_res$p, method = mtc_method)

    # Mark significant rows with asterixes
    ast_q <- rep("", dim(total_res)[1])
    ast_q[total_res$q < 0.05] <- "*"
    total_res$ast_q <- ast_q

    # GENERATE THE PLOT
    total_res$sd_from_mean[total_res$sd_from_mean < 0] <- 0
    graph_theme <- theme_bw(base_size = 12, base_family = "Helvetica") +
        theme(
            panel.grid.major = element_line(size = .5, color = "grey"),
            axis.line = element_line(size = .7, color = "black"), 
                text = element_text(size = 14),
            axis.title.y = element_text(vjust = 0.6)
        ) # + theme(legend.position="none")

    # total_res$
    upperLim <- max(abs(total_res$sd_from_mean))

    total_res$y_ast <- total_res$sd_from_mean * 1.05

    total_res$abs_sd <- abs(total_res$sd_from_mean)

    if ("Direction" %in% colnames(total_res)) {
        the_plot <- ggplot(total_res) +
            geom_bar(aes_string(x = "CellType", y = "abs_sd", 
                                    fill = "Direction"), 
                        position = "dodge", stat = "identity") +
            graph_theme
    } else {
        the_plot <- ggplot(total_res) +
            geom_bar(aes_string(x = "CellType", y = "abs_sd"), 
                        fill = "red", stat = "identity") +
            graph_theme +
            theme(legend.position = "none")
    }

    # Setup the main plot
    the_plot <- the_plot +
        theme(plot.margin = unit(c(1, 0, 0, 0), "mm"), 
                axis.text.x = element_text(angle = 55, hjust = 1)) +
        theme(panel.border = element_rect(colour = "black", 
                                            fill = NA, size = 1)) +
        xlab("") +
        theme(strip.text.y = element_text(angle = 0)) +
        coord_cartesian(ylim = c(0, 1.1 * upperLim)) +
        ylab("Std.Devs. from the mean") + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

    the_plot <- the_plot + 
        scale_y_continuous(breaks = c(0, ceiling(upperLim * 0.66))) + 
        geom_text(aes_string(label = "ast_q", x = "CellType", y = "y_ast"), 
                    size = 10)
    if (multiList) {
        the_plot <- the_plot + 
            facet_grid("list ~ .", scales = "free", space = "free_x")
    }

    # Prepare output
    output <- list()
    output$plain <- the_plot

    if (make_dendro) {
        the_dendrogram <- 
            ctd[[annotLevel]]$plotting$ggdendro_horizontal + 
                theme(plot.margin = unit(c(0, 0, 0, 0), units = "cm"))
        combined_plot <- 
            cowplot::plot_grid(the_dendrogram, the_plot, align = "hv", 
                                ncol = 1, rel_heights = c(1, 1))
        output$withDendro <- combined_plot
    }

    return(output)
}
