#' Check SingleCellExperiment
#'
#' Check whether \code{exp} is a SingleCellExperiment (SCE) object and extract
#' the releveant components.
#'
#' @return List of extracted SCE components.
#'
#' @keywords internal
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays assayNames colData
check_sce <- function(exp,
                      verbose = TRUE) {
    requireNamespace("SummarizedExperiment")
    if (methods::is(exp, "SummarizedExperiment")) {
        # update exp to hold the counts from the SCE
        SE_exp <- exp
        if (!"counts" %in% names(SummarizedExperiment::assays(SE_exp))) {
            if ("raw" %in% names(SummarizedExperiment::assays(SE_exp))) {
                messager("Renaming assay: raw --> counts", v = verbose)
                SummarizedExperiment::assayNames(SE_exp) <-
                    gsub(
                        "raw", "counts",
                        SummarizedExperiment::assayNames(SE_exp)
                    )
            } else {
                stop(
                    "Please ensure counts is the assay name for your raw ",
                    "experiment data in your SE/SCE object"
                )
            }
        }
        exp <- SummarizedExperiment::assays(SE_exp)$counts
        metadata <- SummarizedExperiment::colData(SE_exp)
        # set boolean for later operations
        SE_obj <- TRUE
    } else {
        SE_exp <- NULL
        SE_obj <- FALSE
        metadata <- NULL
    }
    return(list(
        exp = exp,
        metadata = metadata,
        SE_exp = SE_exp,
        SE_obj = SE_obj
    ))
}
