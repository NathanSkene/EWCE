sce_lists_apply <- function(SCE_lists,
                            return_genes = FALSE,
                            level = 2,
                            as_matrix = FALSE,
                            as_DelayedArray = FALSE) {
    lapply(names(SCE_lists), function(x, lvl = level, genes = return_genes) {
        print(x)
        sce_list <- SCE_lists[[x]]
        if (length(sce_list) < lvl) lvl <- length(sce_list)
        sce_lvl <- sce_list[[lvl]]
        print(paste(dim(sce_lvl), collapse = " x "))
        if (as_matrix) {
            for (ass in SummarizedExperiment::assayNames(sce_lvl)) {
                SummarizedExperiment::assay(sce_lvl, ass) <-
                    as(SummarizedExperiment::assay(sce_lvl, ass), "matrix")
            }
        }
        if (as_DelayedArray) {
            for (ass in SummarizedExperiment::assayNames(sce_lvl)) {
                SummarizedExperiment::assay(sce_lvl, ass) <-
                    DelayedArray::DelayedArray(
                        methods::as(SummarizedExperiment::assay(
                            sce_lvl, ass
                        ), "sparseMatrix")
                    )
            }
        }
        if (genes) {
            return(rownames(sce_lvl))
        } else {
            return(sce_lvl)
        }
    }) %>% `names<-`(names(SCE_lists))
}
