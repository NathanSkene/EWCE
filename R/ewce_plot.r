#' Plot EWCE results
#'
#' \code{ewce_plot} generates plots of EWCE enrichment results
#'
#' @param total_res Results data.frame generated using
#' \link[EWCE]{bootstrap_enrichment_test} or
#' \link[EWCE]{ewce_expression_data} functions.
#' Multiple results tables can be
#' merged into one results table, as long as the 'list' column is set to
#' distinguish them.
#' Multiple testing correction is then applied across all merged results.
#' @param mtc_method Method to be used for multiple testing correction.
#' Argument is passed to \link[stats]{p.adjust} (DEFAULT: "bonferroni).
#' @param ctd CellTypeDataset object.
#' Should be provided so that the dendrogram can be taken from it
#' and added to plots
#' @inheritParams cowplot::plot_grid
#'
#' @return A ggplot containing the plot
#'
#' @examples
#' ## Bootstrap significance test,
#' ##  no control for transcript length or GC content
#' ## Use pre-computed results to speed up example
#' full_results <- EWCE::example_bootstrap_results()
#'
#' ## Generate the plot
#' print(EWCE::ewce_plot(
#'     total_res = full_results$results,
#'     mtc_method = "BH"
#' ))
#' @export
#' @import ggplot2
#' @importFrom stats p.adjust
ewce_plot <- function(total_res,
    mtc_method = "bonferroni",
    ctd = NULL,
    align = "v",
    rel_heights = c(.3, 1),
    axis = "lr") {
    requireNamespace("cowplot")
    requireNamespace("gridExtra")
    requireNamespace("grid")
    err_msg <- paste0(
        "ERROR: Invalid mtc_method argument. Please see",
        " '?p.adjust' for valid methods."
    )
    if (!mtc_method %in% c(
        "holm", "hochberg", "hommel",
        "bonferroni", "BH", "BY", "fdr", "none"
    )) {
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
        cells_in_ctd <- function(ctdIN, cells) {
            if (sum(!cells %in% colnames(ctdIN$specificity) == 0)) {
                return(1)
            } else {
                return(0)
            }
        }
        if (length(ctd[[1]]$plotting) > 0) {
            annotLevel <-
                which(unlist(lapply(ctd,
                    FUN = cells_in_ctd,
                    cells = as.character(
                        total_res$CellType
                    )
                )) == 1)
            err_msg2 <- paste0(
                "All of the cells within total_res should come",
                " from a single annotation layer of the CTD"
            )
            if (length(annotLevel) == 0) {
                stop(err_msg2)
            }
        }

        # Set order of cells
        if (length(ctd[[annotLevel]]$plotting) > 0) {
            total_res$CellType <-
                factor(total_res$CellType,
                    levels = ctd[[annotLevel]]$plotting$cell_ordering
                )
        }
    }

    #### Multiple testing correction across all rows ####
    if ("q" %in% colnames(total_res)) {
        total_res$q <- stats::p.adjust(total_res$p,
            method = mtc_method
        )
    }
    #### Mark significant rows with asterixes ####
    ast_q <- rep("", dim(total_res)[1])
    ast_q[total_res$q < 0.05] <- "*"
    total_res$ast_q <- ast_q
    #### Plot ####
    total_res$sd_from_mean[total_res$sd_from_mean < 0] <- 0
    graph_theme <- theme_bw(base_size = 12, base_family = "Helvetica") +
        theme(
            text = element_text(size = 14),
            axis.title.y = element_text(vjust = 0.6),
            strip.background = element_rect(fill = "white"),
            strip.text = element_text(color = "black")
        )

    upperLim <- max(abs(total_res$sd_from_mean), na.rm = TRUE)
    total_res$y_ast <- total_res$sd_from_mean * 1.05
    total_res$abs_sd <- abs(total_res$sd_from_mean)

    if ("Direction" %in% colnames(total_res)) {
        the_plot <- ggplot(total_res) +
            geom_bar(aes_string(
                x = "CellType", y = "abs_sd",
                fill = "Direction"
            ),
            position = "dodge", stat = "identity"
            ) +
            graph_theme
    } else {
        the_plot <- ggplot(total_res) +
            geom_bar(aes_string(x = "CellType", y = "abs_sd", fill = "abs_sd"),
                stat = "identity"
            ) +
            scale_fill_gradient(low = "blue", high = "red") +
            graph_theme +
            theme(legend.position = "none")
    }

    # Setup the main plot
    the_plot <- the_plot +
        theme(
            plot.margin = grid::unit(c(1, 0, 0, 0), "mm"),
            axis.text.x = element_text(angle = 55, hjust = 1)
        ) +
        theme(panel.border = element_rect(
            colour = "black",
            fill = NA, size = 1
        )) +
        xlab("") +
        theme(strip.text.y = element_text(angle = 0)) +
        coord_cartesian(ylim = c(0, 1.1 * upperLim)) +
        ylab("Std.Devs. from the mean") +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

    the_plot <- the_plot +
        scale_y_continuous(breaks = c(0, ceiling(upperLim * 0.66))) +
        geom_text(aes_string(label = "ast_q", x = "CellType", y = "y_ast"),
            size = 10
        )
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
            theme(plot.margin = unit(c(0, 0, 0, 0), units = "cm")) +
            scale_x_discrete(breaks = total_res$CellType)
        # scale_x_discrete to set the mapping of the dendro to the x axis scale
        combined_plot <-
            cowplot::plot_grid(the_dendrogram, the_plot,
                align = align,
                rel_heights = rel_heights,
                axis = axis,
                ncol = 1
            )
        # align arg to "v" and rel_heights c(0.3,1) make dend & barchart closer
        # two prev comments adjustments by Robert Gordon-Smith
        output$withDendro <- combined_plot
    }

    return(output)
}
