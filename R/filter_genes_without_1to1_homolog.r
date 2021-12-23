#' filter_genes_without_1to1_homolog
#'
#' Deprecated function. Please use \link[EWCE]{filter_nonorthologs} instead.
#'
#' @inherit filter_nonorthologs
#' @export
filter_genes_without_1to1_homolog <- function(filenames,
                                              input_species = "mouse",
                                              convert_nonhuman_genes = TRUE,
                                              annot_levels = NULL,
                                              suffix = "_orthologs",
                                              verbose = TRUE) {
    .Deprecated("filter_nonorthologs")
    newFilenames <- filter_nonorthologs(filenames,
        input_species = "mouse",
        convert_nonhuman_genes = TRUE,
        annot_levels = NULL,
        suffix = "_orthologs",
        verbose = TRUE
    )
    return(newFilenames)
}
