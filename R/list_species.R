#' List all species
#'
#' List all species that EWCE can convert genes from/to.
#' Wrapper function for \link[orthogene]{map_species}.
#'
#' @param verbose Print messages.
#'
#' @return List of species EWCE can input/output genes as.
#'
#' @export
#' @importFrom orthogene map_species
#' @examples
#' list_species()
list_species <- function(verbose = TRUE) {
    orthogene::map_species(
        species = NULL,
        method = "homologene",
        verbose = verbose
    )
}
