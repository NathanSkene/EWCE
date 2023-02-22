#' sce_merged_apply
#'
#' Merge a list of SingleCellExperiments.
#'
#' @return Merged SingleCellExperiment.
#'
#' @keywords internal
#' @importFrom methods as
#' @importFrom SummarizedExperiment assayNames assay
sce_merged_apply <- function(SCE_merged,
                             as_sparse = TRUE,
                             as_DelayedArray = FALSE) {
    lapply(names(SCE_merged), function(lvl,
                                       .as_DelayedArray = as_DelayedArray,
                                       .as_sparse = as_sparse) {
        print(lvl)
        sce_lvl <- SCE_merged[[lvl]]
        if (.as_sparse) {
            for (ass in SummarizedExperiment::assayNames(sce_lvl)) {
                SummarizedExperiment::assay(sce_lvl, ass) <-
                    methods::as(
                        SummarizedExperiment::assay(sce_lvl, ass),
                        "sparseMatrix"
                    )
            }
        }
        if (.as_DelayedArray) {
            for (ass in SummarizedExperiment::assayNames(sce_lvl)) {
                SummarizedExperiment::assay(sce_lvl, ass) <-
                    DelayedArray::DelayedArray(
                        methods::as(
                            SummarizedExperiment::assay(sce_lvl, ass),
                            "sparseMatrix"
                        )
                    )
            }
        }
        return(sce_lvl)
    }) |> `names<-`(names(SCE_merged))
}
