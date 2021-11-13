#' Convert a CellTypeDataset into standardized format
#'
#' This function will take a CTD,
#' drop all genes without 1:1 orthologs with the
#' \code{output_species} ("human" by default),
#' convert the remaining genes to gene symbols,
#' assign names to each level,
#' and convert all matrices to sparse matrices and/or \code{DelayedArray}.
#'
#' @param ctd Input CellTypeData.
#' @param dataset CellTypeData. name.
#' @param force_new_quantiles By default, quantile computation is
#' skipped if they have already been computed.
#' Set \code{=TRUE} to override this and generate new quantiles.
#' @param force_standardise If \code{ctd} has already been standardised, whether
#' to rerun standardisation anyway (Default: \code{FALSE}).
#' @param remove_unlabeled_clusters Remove any samples that have
#'  numeric column names.
#' @param numberOfBins Number of non-zero quantile bins.
#' @param keep_annot Keep the column annotation data if provided.
#' @param keep_plots Keep the dendrograms if provided.
#' @param verbose Print messages.
#' Set \code{verbose=2} if you want to print all messages
#'  from internal functions as well.
#' @inheritParams extract_matrix
#' @inheritParams drop_uninformative_genes
#' @inheritParams orthogene::convert_orthologs
#'
#' @return Standardised CellTypeDataset.
#'
#' @examples
#' ctd <- ewceData::ctd()
#' ctd_std <- standardise_ctd(
#'     ctd = ctd,
#'     input_species = "mouse",
#'     dataset = "Zeisel2016"
#' )
#' @export
standardise_ctd <- function(ctd,
                            dataset,
                            input_species = NULL,
                            output_species = "human",
                            non121_strategy = "drop_both_species",
                            force_new_quantiles = TRUE,
                            force_standardise = FALSE,
                            remove_unlabeled_clusters = FALSE,
                            numberOfBins = 40,
                            keep_annot = TRUE,
                            keep_plots = TRUE,
                            as_sparse = TRUE,
                            as_DelayedArray = FALSE,
                            rename_columns = TRUE,
                            make_columns_unique = FALSE,
                            verbose = TRUE) {
    #### Handle verbosity at different levels ####
    if (verbose == 2) {
        verbose <- TRUE
        verbose2 <- TRUE
    } else {
        verbose2 <- FALSE
    }
    if (is_ctd_standardised(ctd = ctd) && (force_standardise == FALSE)) {
        messager("ctd is already standardised. Returning original ctd.\n",
            "Set force_standardise=TRUE to re-standardise.",
            v = verbose
        )
        return(ctd)
    }
    #### Iterate over CTD levels ####
    messager("Standardising CellTypeDataset", v = verbose)
    new_ctd <- lapply(seq(1, length(ctd)), function(lvl) {
        messager("Level:", lvl, v = verbose)
        mean_exp <- extract_matrix(
            ctd = ctd,
            input_species = input_species,
            output_species = output_species,
            dataset = dataset,
            level = lvl,
            metric = "mean_exp",
            non121_strategy = non121_strategy,
            numberOfBins = numberOfBins,
            remove_unlabeled_clusters = remove_unlabeled_clusters,
            force_new_quantiles = force_new_quantiles,
            as_sparse = as_sparse,
            as_DelayedArray = as_DelayedArray,
            rename_columns = rename_columns,
            make_columns_unique = make_columns_unique,
            verbose = verbose2
        )
        spec <- extract_matrix(
            ctd = ctd,
            input_species = input_species,
            output_species = output_species,
            dataset = dataset,
            level = lvl,
            metric = "specificity",
            non121_strategy = non121_strategy,
            numberOfBins = numberOfBins,
            remove_unlabeled_clusters = remove_unlabeled_clusters,
            force_new_quantiles = force_new_quantiles,
            as_sparse = as_sparse,
            as_DelayedArray = as_DelayedArray,
            verbose = verbose2
        )
        specQ <- extract_matrix(
            ctd = ctd,
            input_species = input_species,
            output_species = output_species,
            dataset = dataset,
            level = lvl,
            metric = "specificity_quantiles",
            non121_strategy = non121_strategy,
            numberOfBins = numberOfBins,
            remove_unlabeled_clusters = remove_unlabeled_clusters,
            force_new_quantiles = force_new_quantiles,
            as_sparse = as_sparse,
            as_DelayedArray = as_DelayedArray,
            verbose = verbose2
        )
        #### Add extra items ####
        annot <- if ("annot" %in% names(ctd[[lvl]]) & keep_annot) {
            ctd[[lvl]]$annot
        } else {
            NULL
        }
        plotting <- if ("plotting" %in% names(ctd[[lvl]]) & keep_plots) {
            ctd[[lvl]]$plotting
        } else {
            NULL
        }
        return(list(
            "mean_exp" = mean_exp,
            "specificity" = spec,
            "specificity_quantiles" = specQ,
            "annot" = annot,
            "plotting" = plotting,
            "standardised" = TRUE,
            "gene_species" = output_species
        ))
    })
    #### Name CTD levels ####
    if (is.null(names(ctd))) {
        names(new_ctd) <- paste0("level_", seq(1, length(ctd)))
    } else {
        names(new_ctd) <- names(ctd)
    }
    return(new_ctd)
}
