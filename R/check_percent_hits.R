#' Get percentage of target cell type hits
#'
#' After you run \link[EWCE]{bootstrap_enrichment_test},
#'  check what percentage of significantly enriched
#'  cell types match an expected cell type.
#'
#' @param boot_res \code{bootstrap_enrichment_test} results.
#' @param target_celltype Substring to search to matching
#'  cell types (case-insensitive).
#' @param mtc_method Multiple-testing correction method.
#' @param q_threshold Corrected significance threshold.
#' @param verbose Print messages.
#'
#' @returns Report list.
#'
#' @export
#' @examples
#' ## Bootstrap significance test,
#' ##  no control for transcript length or GC content
#' ## Use pre-computed results to speed up example
#' full_results <- EWCE::example_bootstrap_results()
#'
#' report <- EWCE::check_percent_hits(
#'     boot_res = full_results,
#'     target_celltype = "microglia"
#' )
check_percent_hits <- function(boot_res,
                               target_celltype,
                               mtc_method = "bonferroni",
                               q_threshold = .05,
                               verbose = TRUE) {
    sig_results <- get_sig_results(
        full_results = boot_res,
        mtc_method = mtc_method,
        q_threshold = q_threshold,
        verbose = verbose
    )
    # if(any(ctd_reference=="Zeisel2018")){
    #     z18.terms <- search_Zeisel2018_celltypes(
    # target_celltype=target_celltype,
    # verbose = F)
    #     target_hits <- grep(paste(z18.terms,collapse = "|"),
    #                         unique(sig_results$CellType), value = TRUE)
    # }else {
    target_hits <- grep(target_celltype, sig_results$CellType,
        ignore.case = TRUE, value = TRUE
    )
    z18.terms <- NULL
    # }
    percent_hits <- round(length(unique(target_hits)) /
        length(unique(sig_results$CellType)) * 100, 1)
    msg <- paste(
        paste0(percent_hits, "%"),
        "of hits are of the target cell type."
    )
    messager(msg, v = verbose)
    return(list(
        target_hits = target_hits,
        percent_hits = percent_hits,
        target_celltype = target_celltype,
        z18.terms = z18.terms
    ))
}
