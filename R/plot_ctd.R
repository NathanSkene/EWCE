#' Plot \emph{CellTypeData} metrics
#'
#' Plot \emph{CellTypeData} metrics such as mean_exp, specificity and/or
#' specificity_quantiles.
#'
#' @param ctd CellTypeDataset.
#' @param genes Which genes in \code{ctd} to plot.
#' @param level Annotation level in \code{ctd} to plot.
#' @param metric Which metric in the \code{ctd} to plot:
#' \itemize{
#' \item{"mean_exp"}
#' \item{"specificity"}
#' \item{"specificity_quantiles"}
#' }
#' @param show_plot Whether to print the plot or simply return it.
#'
#' @return ggplot object.
#'
#' @examples
#' ctd <- ewceData::ctd()
#' plt <- EWCE::plot_ctd(ctd, genes = c("Apoe", "Gfap", "Gapdh"))
#' @export
#' @import ggplot2
#' @importFrom  stringr str_to_sentence
#' @importFrom  reshape2 melt
plot_ctd <- function(ctd,
                     genes,
                     level = 1,
                     metric = "specificity",
                     show_plot = TRUE) {
    #### Standardise metric name ####
    if (tolower(metric) %in% c(
        "expr", "exp", "expression",
        "mean_exp", "avgexp"
    )) {
        metric <- "mean_exp"
    }
    metric <- stringr::str_to_sentence(metric)
    mat <- ctd[[level]][[tolower(metric)]]
    genes <- genes[genes %in% rownames(mat)]
    plot_data <- reshape2::melt(mat[genes, ], id.vars = "genes")
    colnames(plot_data) <- c("Gene", "Celltype", metric)

    gp <- ggplot(
        plot_data,
        aes_string(x = "Celltype", y = metric, fill = metric)
    ) +
        scale_fill_gradient(low = "blue", high = "red") +
        geom_bar(stat = "identity") +
        facet_grid(Gene ~ .) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "white"),
            strip.text = element_text(color = "black")
        )

    if (metric == "Specificity") {
        gp <- gp + scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1))
    }
    if (show_plot) print(gp)
    return(gp)
}
