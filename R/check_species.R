check_species <- function(genelistSpecies = NULL,
                          sctSpecies = NULL,
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
    return(list(
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies
    ))
}
