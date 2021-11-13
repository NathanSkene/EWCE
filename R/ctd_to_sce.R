#' CellTypeDataset to SingleCellExperiment
#'
#' Copied from \href{https://github.com/bschilder/scKirby}{scKirby},
#'  which is not yet on CRAN or Bioconductor.
#'
#' @param object CellTypeDataset object.
#' @param as_sparse Store SingleCellExperiment matrices as sparse.
#' @param as_DelayedArray Store SingleCellExperiment matrices as DelayedArray.
#' @param verbose Print messages.
#'
#' @return SingleCellExperiment
#'
#' @examples
#' ctd <- ewceData::ctd()
#' sce <- EWCE::ctd_to_sce(ctd)
#' @export
#' @importFrom dplyr %>%
ctd_to_sce <- function(object,
                       as_sparse = TRUE,
                       as_DelayedArray = FALSE,
                       verbose = TRUE) {
    messager("+ CTD ==> SingleCellExperiment", v = verbose)
    ctd <- object
    #### Name CTD levels ####
    if (is.null(names(ctd))) {
        names(ctd) <- paste0("level_", seq(1, length(ctd)))
    } else {
        names(ctd) <- names(ctd)
    }
    sce_list <- lapply(names(ctd), function(lvl) {
        messager("Converting level: ", lvl, v = verbose)
        ctd_lvl <- ctd[[lvl]]
        #### Use matrices that are present ###
        matrix_list <- list()
        for (mtx_name in c(
            "mean_exp", "median_exp",
            "specificity", "median_specificity", "specificity_quantiles"
        )) {
            if (mtx_name %in% names(ctd_lvl)) {
                mtx <- ctd_lvl[[mtx_name]]
                mtx <- to_sparse_matrix(
                    exp = mtx,
                    as_sparse = as_sparse,
                    verbose = FALSE
                )
                mtx <- to_delayed_array(
                    exp = mtx,
                    as_DelayedArray = as_DelayedArray,
                    verbose = FALSE
                )
                matrix_list[[mtx_name]] <- mtx
            }
        }
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = matrix_list,
            colData = data.frame(colnames(matrix_list[[1]])) %>%
                `colnames<-`(lvl),
            rowData = data.frame(
                gene = row.names(matrix_list[[1]]),
                row.names = row.names(matrix_list[[1]])
            )
        )
        # sce <- check_sce_rownames(sce, verbose = verbose)
    }) %>% `names<-`(names(ctd))
    ## "SCE_list" class messes up other functions that expect class "list"
    # class(sce_list) <- "SCE_list"
    return(sce_list)
}
