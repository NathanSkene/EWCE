delayedarray_normalize <- function(exp,
                                   log_norm = TRUE,
                                   min_max = TRUE,
                                   plot_hists = FALSE) {
    mat <- exp
    core_allocation <- assign_cores(worker_cores = .90)
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
