#' Extract a matrix from a CellTypeDataset
#'
#' Extracts a particular matrix (e.g., mean_exp, specificity)
#' from a CellTypeDataset object.
#'
#' @param level CTD level to extract from.
#' @param metric Name of the matrix to extract.
#' @param as_sparse Convert to sparse matrix.
#' @param as_DelayedArray Convert to \code{DelayedArray}.
#' @param rename_columns Add species and dataset to column names.
#' @inheritParams standardise_ctd
#' @inheritParams drop_uninformative_genes
#' @inheritParams orthogene::convert_orthologs
#'
#' @return (specificity) matrix.
#'
#' @keywords internal
#' @importFrom DelayedArray rowMaxs DelayedArray
#' @importFrom dplyr %>% mutate group_by slice_max
#' @importFrom methods as
extract_matrix <- function(ctd,
                           dataset,
                           level = 1,
                           input_species = NULL,
                           output_species = "human",
                           metric = "specificity",
                           non121_strategy = "drop_both_species",
                           numberOfBins = 40,
                           remove_unlabeled_clusters = FALSE,
                           force_new_quantiles = FALSE,
                           as_sparse = TRUE,
                           as_DelayedArray = FALSE,
                           rename_columns = FALSE,
                           verbose = TRUE) {
    ### Avoid confusing Biocheck
    Gene <- max_exp <- NULL

    messager("Extracting ", metric, v = verbose)
    ##### Check quantiles ####
    # Check if specificity quantiles are present, and if not compute them.
    if ((metric == "specificity_quantiles") &
        (!"specificity_quantiles" %in% names(ctd[[1]]) |
            force_new_quantiles)) {
        ctd[[level]] <- bin_specificity_into_quantiles(
            ctdIN = ctd[[level]],
            numberOfBins = numberOfBins
        )
    }
    specificity <- ctd[[level]][[metric]]
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
                        colnames(specificity)
                    )
                )
            )
        specificity <- specificity[, labeled_cols]
    }
    #### Ensure SYMBOL format ####
    ## If most of the genes start with an ensembl-style name...
    if (sum(startsWith(rownames(specificity)[seq(1, 100)], "ENS")) > 50) {
        messager("ENSEMBL IDs detected. Converting to HGNC.",
            v = verbose
        )
        specificity <- orthogene::aggregate_mapped_genes(
            gene_df = specificity,
            input_species = input_species,
            non121_strategy = "sum",
            verbose = verbose
        )
    }
    #### Convert orthologs ####
    if (input_species != output_species) {
        specificity <- orthogene::convert_orthologs(
            gene_df = specificity,
            gene_input = "rownames",
            input_species = input_species,
            output_species = output_species,
            non121_strategy = non121_strategy,
            method = "homologene",
            verbose = FALSE
        )
    }
    #### Add input_species/dataset info to column names ####
    if (rename_columns) {
        #### Standardize colnames ####
        colnames(specificity) <- gsub("[.]", "_", colnames(specificity))
        col_prefix <- paste(c(input_species, dataset), collapse = ".")
        colnames(specificity) <- gsub(
            " ", "_",
            paste(col_prefix,
                colnames(specificity),
                sep = "."
            )
        )
    }
    #### Convert to sparse matrix ####
    exp <- to_sparse_matrix(
        exp = specificity,
        as_sparse = as_sparse,
        verbose = verbose
    )
    #### Convert to DelayedArray ####
    specificity <- to_delayed_array(
        exp = specificity,
        as_DelayedArray = as_DelayedArray,
        verbose = verbose
    )
    #### Report ####
    messager(
        "Matrix dimensions:",
        paste(dim(specificity), collapse = " x "),
        v = verbose
    )
    return(specificity)
}
