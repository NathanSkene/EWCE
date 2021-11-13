#' Get summed proportions
#'
#' \code{get_summed_proportions} Given the target gene set, randomly sample
#' gene lists of equal length, obtain the specificity of these and then
#' obtain the mean specificity in each sampled list (and the target list).
#'
#' See \link[EWCE]{bootstrap_enrichment_test} for examples.
#'
#' @param hitGenes list of gene names. The target gene set.
#' @param control_network If \code{geneSizeControl=TRUE},
#' then must provide the control network.
#' @inheritParams bootstrap_enrichment_test
#'
#' @returns A list containing three data frames:
#' \itemize{
#'   \item \code{hit.cells}: vector containing the summed proportion of
#'   expression in each cell type for the target list
#'   \item \code{bootstrap_data}: matrix in which each row represents the
#'   summed proportion of expression in each cell type for one of the
#'   random lists
#'   \item \code{controlledCT}: the controlled cell type (if applicable)
#' }
#'
#' @keywords internal
#' @importFrom data.table data.table rbindlist
#' @importFrom dplyr %>%
#' @importFrom parallel mclapply
get_summed_proportions <- function(hitGenes,
                                   sct_data,
                                   annotLevel,
                                   reps,
                                   no_cores = 1,
                                   geneSizeControl,
                                   controlledCT = NULL,
                                   control_network = NULL,
                                   verbose = TRUE) {
    messager("Computing summed proportions.", v = verbose)
    controlledCT <- fix_celltype_names(celltypes = controlledCT)
    combinedGenes <- rownames(sct_data[[annotLevel]]$mean_exp)
    hitGenes <- hitGenes[hitGenes %in% combinedGenes]
    hit.cells <- cell_list_dist(
        hitGenes = hitGenes,
        sct_data = sct_data,
        annotLevel = annotLevel
    )
    # cell_list_dist gets the summed proportion of 'hitGenes' across all
    # cell types at annotLevel

    # Check control_network provided if geneSizeControl=TRUE
    err_msg <- paste0(
        "ERROR: if geneSizeControl==TRUE then",
        " get_summed_proportions must be passed",
        " control_network as an argument"
    )
    if (geneSizeControl == TRUE) {
        if (is.null(control_network)) {
            stop(err_msg)
        }
    }
    # If controlling for expression of another cell type, then the samples of
    # bootstrapped genes must match the deciles of expression specificity
    if (!is.null(controlledCT)) {
        controlled_bootstrap_set <-
            generate_controlled_bootstrap_geneset(
                hitGenes = hitGenes,
                sct_data = sct_data,
                annotLevel = annotLevel,
                reps = reps,
                controlledCT = controlledCT
            )
    }
    #### Parallelise bootstrapping ####
    bootstrap_data <- parallel::mclapply(seq_len(reps), function(s) {
        # Get 'bootstrap_set'...a list of genes of equivalent length as hitGenes
        if (isTRUE(geneSizeControl)) {
            bootstrap_set <- control_network[s, ]
        } else {
            if (is.null(controlledCT)) {
                bootstrap_set <- sample(
                    combinedGenes,
                    length(hitGenes)
                )
            } else {
                bootstrap_set <- controlled_bootstrap_set[, s]
            }
        }
        # 'bootstrap_data' is a matrix of the summed proportions
        bootstrap_res <- cell_list_dist(
            hitGenes = bootstrap_set,
            sct_data = sct_data,
            annotLevel = annotLevel
        )
        return(data.table::data.table(t(bootstrap_res)))
    }, mc.cores = no_cores) %>%
        data.table::rbindlist()
    #### Conver to matrix format ####
    bootstrap_data <- as.matrix(bootstrap_data)
    return(list(
        hit.cells = hit.cells,
        bootstrap_data = bootstrap_data,
        controlledCT = controlledCT
    ))
}
