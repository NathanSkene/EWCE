#' Example bootstrap celltype enrichment test for transcriptome data
#'
#' Example celltype enrichment
#' results produced by \link[EWCE]{ewce_expression_data}.
#'
#' @param verbose Print messages.
#' @param localHub If working offline, add argument localHub=TRUE to work 
#' with a local, non-updated hub; It will only have resources available that
#' have previously been downloaded. If offline, Please also see BiocManager
#' vignette section on offline use to ensure proper functionality. 
#'
#' @source
#' ## Load the single cell data
#'
#' ctd <- ewceData::ctd()
#'
#' ## Set the parameters for the analysis
#'
#' ## Use 3 bootstrap lists for speed, for publishable analysis use >10,000
#'
#' reps <- 3
#'
#' annotLevel <- 1 # <- Use cell level annotations (i.e. Interneurons)
#'
#' ## Use 5 up/down regulated genes (thresh) for speed, default is 250
#'
#' thresh <- 5
#'
#' ## Load the top table
#'
#' tt_alzh <- ewceData::tt_alzh()
#'
#' tt_results <- EWCE::ewce_expression_data(
#'     sct_data = ctd,
#'     tt = tt_alzh,
#'     annotLevel = 1,
#'     thresh = thresh,
#'     reps = reps,
#'     ttSpecies = "human",
#'     sctSpecies = "mouse"
#' )
#'
#' save(tt_results, file = "inst/extdata/tt_results.rda")
#' @returns List with 5 items.
#' @export
#' @examples
#' tt_results <- EWCE::example_transcriptome_results()
example_transcriptome_results <- function(verbose = TRUE, localHub=FALSE) {
    fname <- system.file("extdata/tt_results.rda",
                         package = "EWCE"
    )
    if (file.exists(fname)) {
        messager("Loading precomputed example transcriptome results.",
                 v = verbose)
        tt_results <- load_rdata(fname)
    } else {
        messager("Recomputing example transcriptome results.",
                 v = verbose)
        ctd <- ewceData::ctd(localHub = localHub)
        reps <- 3
        annotLevel <- 1 # <- Use cell level annotations (i.e. Interneurons)
        thresh <- 5
        ## Load the top table
        tt_alzh <- ewceData::tt_alzh(localHub = localHub)
        tt_results <- ewce_expression_data(
            sct_data = ctd,
            tt = tt_alzh,
            annotLevel = 1,
            thresh = thresh,
            reps = reps,
            ttSpecies = "human",
            sctSpecies = "mouse",
            localHub = localHub
        )
    }
    return(tt_results)
}
