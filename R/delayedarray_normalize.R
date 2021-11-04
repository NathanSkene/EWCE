#' Efficiently normalize a DelayedArray
#'
#' The following is a matrix normalization procedure that takes advantage of
#' functions designed to be more efficient for DelayedArray objects.
#'
#' @param exp Input matrix (e.g. gene expression).
#' @param log_norm Whether to first log-normalise \code{exp}
#'  with \link[base]{log1p}.
#' @param min_max Whether to min/max-normalise \code{exp}.
#' @param no_cores Number of cores to parallelise across.
#'
#' @return Normalised matrix.
#'
#' @keywords internal
delayedarray_normalize <- function(exp,
                                   log_norm = TRUE,
                                   min_max = TRUE,
                                   plot_hists = FALSE,
                                   no_cores = 1) {
    requireNamespace("DelayedArray")
    requireNamespace("graphics")
    mat <- exp
    core_allocation <- assign_cores(worker_cores = no_cores)
    if (log_norm) {
        mat_log <- log1p(mat)
        mat <- mat_log
    }
    if (min_max) {
        col_max <- DelayedArray::colMaxs(mat, na.rm = TRUE)
        col_min <- DelayedArray::colMins(mat, na.rm = TRUE)
        mat_normed <- DelayedArray::t((DelayedArray::t(mat) - col_min) /
            (col_max - col_min))
        mat <- mat_normed
    }
    if (plot_hists) {
        graphics::hist(DelayedArray::colMeans(exp, na.rm = TRUE)) %>% print()
    }
    return(mat)
}
