#' Construct SingleCellExperiment file path
#'
#' Construct a save path for a SingleCellExperiment (SCE).
#'
#' @return SingleCellExperiment file path.
#'
#' @keywords internal
#' @importFrom DelayedArray seed
#' @importFrom SummarizedExperiment assay
#' @importFrom methods slotNames
sce_filepath <- function(sce,
    sce_save_dir = NULL) {
    file_info <- DelayedArray::seed(SummarizedExperiment::assay(sce))
    if ("filepath" %in% methods::slotNames(file_info)) {
        return(file_info@filepath)
    } else {
        return(NULL)
    }
}
