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
#' and added to plots.
#' @param annotLevel An integer indicating which level of \code{ctd} to
#' analyse (\emph{Default: 1}).
#' @param make_dendro Add a dendrogram (requires \code{ctd}).
#' @param heights The relative heights row in the grid. 
#' Will get repeated to match the dimensions of the grid.
#' Passed to \link[patchwork]{wrap_plots}.
#' @inheritParams check_percent_hits 
#'
#' @returns A named list containing versions of the \link[ggplot2]{ggplot}
#'  with and without the dendrogram. 
#'
#' @export
#' @import ggplot2
#' @importFrom stats p.adjust
#' @examples
#' ## Bootstrap significance test,
#' ##  no control for transcript length or GC content
#' ## Use pre-computed results to speed up example
#' total_res <- EWCE::example_bootstrap_results()$results 
#' plt <- ewce_plot(total_res = total_res)
ewce_plot <- function(total_res,
                      mtc_method = "bonferroni",
                      q_threshold = 0.05,
                      ctd = NULL,
                      annotLevel = 1, 
                      heights = c(.3, 1), 
                      make_dendro = FALSE,
                      verbose = TRUE) {
    # templateR:::args2vars(ewce_plot)
     
    requireNamespace("ggplot2")
    requireNamespace("patchwork")    
    
    check_mtc_method(mtc_method = mtc_method)
    multiList <- TRUE
    if (is.null(total_res$list)) multiList <- FALSE
    #### If using dendrogram ####
    if(isTRUE(make_dendro)){
        #### Check if ctd is provided ####
        if(is.null(ctd)){
            messager(
                "Warning: Can only add the dendrogram when ctd is provided.",
                "Setting make_dendro=FALSE.",
                v=verbose)
            make_dendro <- FALSE
        } else {
            # Find the relevant level of the CTD annotation 
            if (length(ctd[[annotLevel]]$plotting) > 0) {
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
            #### Set order of cells ####
            if (length(ctd[[annotLevel]]$plotting) > 0) {
                total_res$CellType <-
                    factor(x = fix_celltype_names(total_res$CellType),
                           levels = fix_celltype_names(
                               ctd[[annotLevel]]$plotting$cell_ordering
                           ), 
                           ordered = TRUE
                    )
            }
        } 
    }
    #### Multiple testing correction across all rows ####
    if (!"q" %in% colnames(total_res)) {
        total_res$q <- stats::p.adjust(total_res$p,
            method = mtc_method
        )
    }
    #### Mark significant rows with asterixes ####
    ast_q <- rep("", dim(total_res)[1])
    ast_q[total_res$q < q_threshold] <- "*"
    total_res$ast_q <- ast_q
    #### Plot ####
    total_res$sd_from_mean[total_res$sd_from_mean < 0] <- 0
    graph_theme <- ggplot2::theme_bw(base_size = 12, 
                                     base_family = "Helvetica") +
        ggplot2::theme(
            text = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(vjust = 0.6),
            strip.background = ggplot2::element_rect(fill = "white"),
            strip.text = ggplot2::element_text(color = "black")
        )

    upperLim <- max(abs(total_res$sd_from_mean), na.rm = TRUE)
    total_res$y_ast <- total_res$sd_from_mean * 1.05
    total_res$abs_sd <- abs(total_res$sd_from_mean)

    if ("Direction" %in% colnames(total_res)) {
        the_plot <- ggplot2::ggplot(total_res) +
            ggplot2::geom_bar(
                ggplot2::aes_string(x = "CellType", y = "abs_sd",
                                    fill = "Direction"
            ),
            position = "dodge", stat = "identity"
            ) +
            graph_theme
    } else {
        the_plot <- ggplot2::ggplot(total_res) +
            ggplot2::geom_bar(
                ggplot2::aes_string(x = "CellType", y = "abs_sd", 
                                    fill = "abs_sd"),
                stat = "identity"
            ) +
            ggplot2::scale_fill_gradient(low = "blue", high = "red") +
            graph_theme +
            ggplot2::theme(legend.position = "none")
    }

    # Setup the main plot
    the_plot <- the_plot +
        ggplot2::theme(
            plot.margin = ggplot2::unit(c(.5, 0, 0, 0), "mm"),
            axis.text.x = ggplot2::element_text(angle = 55, hjust = 1)
        ) +
        ggplot2::theme(panel.border = ggplot2::element_rect(
            colour = "black",
            fill = NA, linewidth = 1
        )) +
        ggplot2::xlab("Cell type") +
        ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) +
        ggplot2::ylab("Std.Devs. from the mean") 

    the_plot <- the_plot +
        ggplot2::scale_y_continuous(breaks = c(0, ceiling(upperLim * 0.66)),
                                    expand = c(0, 1.1)) +
        ggplot2::geom_text(
            ggplot2::aes_string(label = "ast_q", x = "CellType", y = "y_ast"),
            size = 10
        )
    if (isTRUE(multiList)) {
        the_plot <- the_plot +
            ggplot2::facet_grid("list ~ .", 
                                scales = "free", 
                                space = "free_x")
    }
    #### Prepare output list ####
    output <- list()
    output$plain <- the_plot 
    if (isTRUE(make_dendro)) {
        ctdIN <- prep_dendro(ctdIN = ctd[[annotLevel]], 
                             expand = c(0, .66))  
        output$withDendro <- patchwork::wrap_plots(
            ctdIN$plotting$ggdendro_horizontal,
            the_plot, 
            heights = heights,
            ncol = 1)
    }

    return(output)
}
