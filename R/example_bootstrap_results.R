#' Example bootstrap enrichment results
#'
#' Example cell type enrichment
#' results produced by \link[EWCE]{bootstrap_enrichment_test}.
#'
#' @param verbose Print messages.
#'
#' @source
#' \code{
#' # Load the single cell data
#' ctd <- ewceData::ctd()
#' # Set the parameters for the analysis
#' # Use 3 bootstrap lists for speed, for publishable analysis use >=10,000
#' reps <- 3
#' # Load gene list from Alzheimer's disease GWAS
#' example_genelist <- ewceData::example_genelist()
#'
#' # Bootstrap significance test, no control for transcript length or GC content
#' full_results <- EWCE::bootstrap_enrichment_test(
#'     sct_data = ctd,
#'     hits = example_genelist,
#'     reps = reps,
#'     annotLevel = 1,
#'     sctSpecies = "mouse",
#'     genelistSpecies = "human"
#' )
#' bootstrap_results <- full_results
#' save(bootstrap_results,file = "inst/extdata/bootstrap_results.rda")
#' }
#' @returns List with 3 items.
#' @export
#' @examples
#' full_results <- EWCE::example_bootstrap_results()
example_bootstrap_results <- function(verbose = TRUE) {
    fname <- system.file("extdata/bootstrap_results.rda",
        package = "EWCE"
    )
    if (file.exists(fname)) {
        messager("Loading precomputed example bootstrap results.", v = verbose)
        full_results <- load_rdata(fname)
    } else {
        messager("Recomputing example bootstrap results.", v = verbose)
        ctd <- ewceData::ctd()
        hits <- ewceData::example_genelist()
        full_results <- bootstrap_enrichment_test(
            sct_data = ctd,
            hits = hits,
            reps = 100,
            annotLevel = 1,
            sctSpecies = "mouse",
            genelistSpecies = "human",
            verbose = verbose
        )
    }
    return(full_results)
}
