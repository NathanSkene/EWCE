#' Extract a matrix from a CellTypeDataset
#'
#' Extracts a particular matrix (e.g., mean_exp, specificity)
#' from a CellTypeDataset object.
#'
#' @param level CTD level to extract from.
#' @param metric Name of the matrix to extract.
#' @param as_sparse Convert to sparse matrix.
#' @param as_DelayedArray Convert to \code{DelayedArray}.
#' @param rename_columns Remove \code{replace_chars} from column names.
#' @param make_columns_unique Rename each columns with the prefix
#'  \code{dataset.species.celltype}.
#' @inheritParams standardise_ctd
#' @inheritParams drop_uninformative_genes
#' @inheritParams orthogene::convert_orthologs
#' @inheritDotParams orthogene::convert_orthologs
#'
#' @return (specificity) matrix.
#'
#' @keywords internal
#' @importFrom DelayedArray rowMaxs DelayedArray
#' @importFrom dplyr mutate group_by slice_max
#' @importFrom methods as
extract_matrix <- function(ctd,
                           dataset,
                           level = 1,
                           input_species = NULL,
                           output_species = "human",
                           metric = "specificity",
                           non121_strategy = "drop_both_species",
                           method = "homologene",
                           numberOfBins = 40,
                           remove_unlabeled_clusters = FALSE,
                           force_new_quantiles = FALSE,
                           as_sparse = TRUE,
                           as_DelayedArray = FALSE,
                           rename_columns = TRUE,
                           make_columns_unique = FALSE,
                           verbose = TRUE,
                           ...) {
    ### Avoid confusing Biocheck
    Gene <- max_exp <- NULL

    messager("Extracting ", metric, v = verbose)
    ##### Check quantiles ####
    # Check if specificity quantiles are present, and if not compute them.
    if ((metric == "specificity_quantiles") &&
        (!"specificity_quantiles" %in% names(ctd[[1]]) |
            force_new_quantiles)) {
        ctd[[level]] <- bin_specificity_into_quantiles(
            ctdIN = ctd[[level]],
            numberOfBins = numberOfBins
        )
    }
    mat <- ctd[[level]][[metric]]
    #### Remove unlabeled clusters ####
    if (remove_unlabeled_clusters) {
        messager("+ Removing unlabeled cell types.", v = verbose)
        labeled_cols <-
            is.na(
                as.numeric(
                    gsub(
                        paste(c(input_species, dataset, "_"),
                            collapse = "|"
                        ), "",
                        colnames(mat)
                    )
                )
            )
        mat <- mat[, labeled_cols]
    }
    #### Ensure SYMBOL format ####
    ## If most of the genes start with an ensembl-style name...
    if (sum(startsWith(rownames(mat)[seq(1, 100)], "ENS")) > 50) {
        messager("ENSEMBL IDs detected. Converting to HGNC.",
            v = verbose
        )
        mat <- orthogene::aggregate_mapped_genes(
            gene_df = mat,
            input_species = input_species, 
            verbose = verbose
        )
    }
    #### Convert orthologs ####
    if (input_species != output_species) {
        mat <- orthogene::convert_orthologs(
            gene_df = mat,
            gene_input = "rownames",
            input_species = input_species,
            output_species = output_species,
            non121_strategy = non121_strategy,
            method = method,
            verbose = FALSE,
            ...
        )
    }
    #### Replace problematic characters in colnames ####
    if (rename_columns) {
        colnames(mat) <- fix_celltype_names(celltypes = colnames(mat))
    }
    #### Make each colname unique (when merging multiple CTD) ####
    if (make_columns_unique) {
        col_prefix <- paste(c(input_species, dataset), collapse = ".")
        colnames(mat) <- make.unique(
            paste(col_prefix, colnames(mat), sep = ".")
        )
    }
    #### Convert to sparse matrix ####
    mat <- to_sparse_matrix(
        exp = mat,
        as_sparse = as_sparse,
        verbose = verbose
    )
    #### Convert to DelayedArray ####
    mat <- to_delayed_array(
        exp = mat,
        as_DelayedArray = as_DelayedArray,
        verbose = verbose
    )
    #### Report ####
    messager(
        "Matrix dimensions:",
        paste(dim(mat), collapse = " x "),
        v = verbose
    )
    return(mat)
}
