#' Run DGE:  \pkg{MAST}
#'
#' Run Differential Gene Expression with \pkg{MAST}.
#'
#' @param no_cores Number of cores to parallelise DGE across.
#' @inheritParams drop_uninformative_genes
#' @inheritParams MAST::zlm
#'
#' @source \href{https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html}{MAST tutorial}
#'
#' @return \code{MAST} results
#'
#' @keywords internal
run_mast <- function(exp,
                     level2annot,
                     test = "LRT",
                     mtc_method = "BH",
                     no_cores = 1,
                     ...) {
    requireNamespace("MAST")
    sca <- MAST::FromMatrix(
        exprsArray = exp,
        cData = data.frame(Class = level2annot),
        fData = data.frame(Gene = rownames(exp)),
        check_sanity = FALSE
    )
    options(mc.cores = no_cores)
    if (test == "LRT") {
        mast_res <- MAST::LRT(
            sca = sca,
            comparison = "Class",
            ...
        )
    } # else {
    #     zlm_res <- MAST::zlm(formula = ~Class,
    #                      sca = sca,
    #                      parallel = no_cores > 1,
    #                      ...
    #                      )
    #     summ <- summary(object = zlm_res,
    #                     doLRT = "Class")
    #     summaryDt <- summ$datatable
    #     fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
    #                       summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    #
    #     fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    #     fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
    #     setorder(fcHurdleSig, fdr)
    #
    #     # mast_res <- MAST::waldTest(object = r,
    #     #                            hypothesis = MAST::CoefficientHypothesis("Class"))
    # }
    mast_res$q <- stats::p.adjust(
        p = mast_res$p.value,
        method = mtc_method
    )
    return(mast_res)
}
