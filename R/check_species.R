#' Check species 
#' 
#' If species arguments are \code{NULL}, set default species.
#' 
#' @inheritParams bootstrap_enrichment_test
#' @returns List of corrected species names.
#' 
#' @keywords internal
check_species <- function(genelistSpecies = NULL,
                          sctSpecies = NULL,
                          sctSpecies_origin = NULL,
                          verbose = TRUE) {
    if (is.null(genelistSpecies)) {
        messager(
            "Warning: genelistSpecies not provided.",
            "Setting to 'human' by default.",
            v = verbose
        )
        genelistSpecies <- "human"
    }
    if (is.null(sctSpecies)) {
        messager(
            "Warning: sctSpecies not provided.",
            "Setting to 'mouse' by default.",
            v = verbose
        )
        sctSpecies <- "mouse"
    }
    if (is.null(sctSpecies_origin)) {
        messager(
            "Warning: sctSpecies_origin not provided.",
            "Setting to 'mouse' by default.",
            v = verbose
        )
        sctSpecies_origin <- "mouse"
    }
    return(list(
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies,
        sctSpecies_origin = sctSpecies_origin
    ))
}
