#' Example bootstrap enrichment results
#'
#' Example cell type enrichment
#' results produced by \link[EWCE]{bootstrap_enrichment_test}.
#'
#' @param verbose Print messages.
#' @param localHub If working offline, add argument localHub=TRUE to work 
#' with a local, non-updated hub; It will only have resources available that
#' have previously been downloaded. If offline, Please also see BiocManager
#' vignette section on offline use to ensure proper functionality. 
#' @source
#' # Load the single cell data
#'
#' ctd <- ewceData::ctd()
#'
#' # Set the parameters for the analysis
#'
#' # Use 3 bootstrap lists for speed, for publishable analysis use >=10,000
#'
#' reps <- 3
#'
#' # Load gene list from Alzheimer's disease GWAS
#'
#' example_genelist <- ewceData::example_genelist()
#'
#' # Bootstrap significance test, no control for transcript length or GC content
#'
#' full_results <- EWCE::bootstrap_enrichment_test(
#'     sct_data = ctd,
#'     hits = example_genelist,
#'     reps = reps,
#'     annotLevel = 1,
#'     sctSpecies = "mouse",
#'     genelistSpecies = "human"
#' )
#'
#' bootstrap_results <- full_results
#'
#' save(bootstrap_results,file = "inst/extdata/bootstrap_results.rda")
#' @returns List with 3 items.
#' @export
#' @examples
#' full_results <- EWCE::example_bootstrap_results()
example_bootstrap_results <- function(verbose = TRUE,localHub = FALSE) {
    fname <- system.file("extdata/bootstrap_results.rda",
        package = "EWCE"
    )
    if (file.exists(fname)) {
        messager("Loading precomputed example bootstrap results.", v = verbose)
        full_results <- load_rdata(fname)
    } else {
        messager("Recomputing example bootstrap results.", v = verbose)
        ctd <- ewceData::ctd(localHub = localHub)
        hits <- ewceData::example_genelist(localHub = localHub)
        full_results <- bootstrap_enrichment_test(
            sct_data = ctd,
            hits = hits,
            reps = 100,
            annotLevel = 1,
            sctSpecies = "mouse",
            genelistSpecies = "human",
            verbose = verbose,
            localHub = localHub
        )
    }
    return(full_results)
}
