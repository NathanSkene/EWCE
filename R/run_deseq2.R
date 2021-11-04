#' Run DGE:  \pkg{DESeq2}
#'
#' Run Differential Gene Expression with \pkg{DESeq2}.
#'
#' @inheritParams drop_uninformative_genes
#' @inheritParams DESeq2::DESeq
#'
#' @return \code{DESeq} results
#'
#' @keywords internal
#' @importFrom stats formula
run_deseq2 <- function(exp,
                        level2annot,
                        test = "LRT",
                        no_cores = 1,
                        verbose = TRUE,
                        ...) {
    requireNamespace("DESeq2")
    messager("DGE:: DESeq2...", v = verbose)
    core_allocation <- assign_cores(worker_cores = no_cores)
    # NOTE:: When you're running DESeq2 on sparse exp data,
    ## there are two ways to avoid issues when DESeq() tries to log your data:
    ## 1) add 1 to you expression matrix (much faster).
    ## 2) set sfType = "iterate" to
    ##    enable iterative size factor estimation (veerrrry slow).
    dds <- DESeq2::DESeqDataSetFromMatrix(exp + 1,
        colData = data.frame(level2annot = level2annot),
        design = stats::formula(paste("~", "level2annot"))
    )
    dds <- DESeq2::DESeq(dds,
        # Best for scRNAseq data.
        test = test,
        reduced = ~1,
        # DESeq2 v1.31.10 (not yet released on BioC)
        # now has glmGamPoi integrated directly!
        ## https://gi th ub.com/mikelove/DESeq2/issues/29
        ## default="parametric"
        # fitType="glmGamPoi",
        # sfType = "iterate",
        parallel = no_cores > 1,
        ...
    )
    dds_res <- DESeq2::results(dds)
    return(dds_res)
}
