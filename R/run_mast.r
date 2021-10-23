#' Run DGE:  \pkg{MAST}
#'
#' Run Differential Gene Expression with \pkg{MAST}.
#'
#' @return \code{MAST} results
#' @inheritParams drop_uninformative_genes
#' @inheritParams MAST::zlm
#' @keywords internal
run_mast <- function(exp,
                     level2annot,
                     no_cores = 1,
                     ...) {
    requireNamespace("MAST")
    sca_delay <- MAST::FromMatrix(
        exprsArray = list(exp),
        cData = level2annot,
        fData = rownames(exp),
        check_sanity = FALSE
    )
    options(mc.cores = no_cores)
    res <- MAST::zlm(~Class,
        sca = sca_delay,
        parallel = no_cores > 1,
        ...
    )
    mast_res <- MAST::waldTest(res, MAST::CoefficientHypothesis("Stim.ConditionUnstim"))
    return(mast_res)
}
