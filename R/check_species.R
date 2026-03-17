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
    validate_species <- function(species, arg_name) {
        species_map <- orthogene::map_species(
            species = species,
            method = "homologene",
            verbose = FALSE
        )
        if (length(species_map) == 0 || all(is.na(species_map))) {
            stop(arg_name, " is not a recognised species: '", species, "'.")
        }
        names(species_map)[1]
    }
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
    genelistSpecies <- validate_species(
        species = genelistSpecies,
        arg_name = "genelistSpecies"
    )
    sctSpecies <- validate_species(
        species = sctSpecies,
        arg_name = "sctSpecies"
    )
    sctSpecies_origin <- validate_species(
        species = sctSpecies_origin,
        arg_name = "sctSpecies_origin"
    )
    return(list(
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies,
        sctSpecies_origin = sctSpecies_origin
    ))
}
