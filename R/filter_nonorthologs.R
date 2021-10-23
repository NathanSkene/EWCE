#' Filter non-orthologs
#'
#' \code{filter_nonorthologs} Takes the filenames of CellTypeData files,
#' loads them,  drops any genes which don't have a 1:1 orthologs with humans,
#' and then convert the gene to human orthologs.
#' The new files are then saved to disk, appending
#' '_orthologs' to the file name.
#'
#' \bold{Note:} This function replaces the original
#'  \code{filter_genes_without_1to1_homolog} function.
#' \code{filter_genes_without_1to1_homolog} is
#' now a wrapper for \code{filter_nonorthologs}.
#'
#' @param filenames List of file names for sct_data saved as \emph{.rda} files.
#' @param input_species Which species the gene names in \code{exp} come from.
#' @param annot_levels [Optional] Names of each annotation level.
#' @param suffix Suffix to add to the file name (right before \emph{.rda}).
#' @param convert_nonhuman_genes Whether to convert the \code{exp}
#' row names to human gene names.
#' @param verbose Print messages.
#'
#' @return List of the filtered CellTypeData file names.
#'
#' @examples
#' # Load the single cell data
#' ctd <- ewceData::ctd()
#' tmp <- tempfile()
#' save(ctd, file = tmp)
#' fNames_ALLCELLS_orths <- EWCE::filter_nonorthologs(filenames = tmp)
#' @export
#' @importFrom orthogene convert_orthologs
filter_nonorthologs <- function(filenames,
                                input_species = NULL,
                                convert_nonhuman_genes = TRUE,
                                annot_levels = NULL,
                                suffix = "_orthologs",
                                verbose = TRUE) {
    if (is.null(input_species)) {
        messager("No input_species provided. Setting to 'mouse'", v = verbose)
        input_species <- "mouse"
    }
    is_filename <- methods::is(filenames[1], "character")
    is_ctd_list <- all(
        c("annot", "mean_exp", "specificity") %in% names(filenames[[1]][[1]])
    )
    if (!is_ctd_list) {
        is_ctd <- is_celltypedataset(ctd = filenames)
        if (is_ctd) {
            filenames <- list(ctd1 = filenames)
            is_ctd_list <- TRUE
        }
    }
    if (is_ctd_list) {
        newFilenames <- list()
    } else {
        newFilenames <- c()
    }
    for (i in seq(1, length(filenames))) {
        #### Read in file ####
        if (is_ctd_list) {
            ctd <- filenames[[i]]
        } else {
            ff <- filenames[[i]]
            ctd <- load_rdata(ff)
        }
        if (is.null(annot_levels)) annot_levels <- seq(1, length(ctd))
        for (lvl in annot_levels) {
            messager("+ Processing level", lvl, "...", v = verbose)
            ctd_slots <- names(ctd[[lvl]])
            ctd_slots <- ctd_slots[!ctd_slots %in% c("annot", "plotting")]
            for (x in ctd_slots) {
                messager("Processing", x, v = verbose)
                ### Convert to human orthologs
                ctd[[lvl]][[x]] <- orthogene::convert_orthologs(
                    gene_df = ctd[[lvl]][[x]],
                    gene_input = "rownames",
                    input_species = input_species,
                    non121_strategy = "drop_both_species",
                    method = "homologene",
                    verbose = FALSE
                )
            }
        }
        if (is_ctd_list) {
            newFilenames[[paste0("ctd", i)]] <- ctd
        } else {
            ff2 <- gsub("\\.rda", paste0(suffix, "\\.rda"), ff)
            save(ctd, file = ff2)
            newFilenames <- c(newFilenames, ff2)
        }
    }
    return(newFilenames)
}
