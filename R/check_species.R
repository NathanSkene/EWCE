#' Check species 
#' 
#' If species arguments are \code{NULL}, set default species.
#' 
#' @param sctSpecies_origin_default Default value for \code{sctSpecies_origin}.
#' @inheritParams bootstrap_enrichment_test
#' @returns List of corrected species names.
#' 
#' @keywords internal
check_species <- function(genelistSpecies = NULL,
                          sctSpecies = NULL,
                          sctSpecies_origin = NULL,
                          sctSpecies_origin_default="mouse",
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
    if (is.null(sctSpecies_origin) &&
        !is.null(sctSpecies_origin_default)) {
        messager(
            "Warning: sctSpecies_origin not provided.",
            "Setting to",shQuote(sctSpecies_origin_default),
            "by default.",
            v = verbose
        )
        sctSpecies_origin <- sctSpecies_origin_default
    }
    return(list(
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies,
        sctSpecies_origin = sctSpecies_origin
    ))
}
