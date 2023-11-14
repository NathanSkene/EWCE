#' Get summed proportions
#'
#' \code{get_summed_proportions} Given the target gene set, randomly sample
#' gene lists of equal length, obtain the specificity of these and then
#' obtain the mean specificity in each sampled list (and the target list).
#'
#' See \link[EWCE]{bootstrap_enrichment_test} for examples.
#'
#' @param hits list of gene names. The target gene set.
#' @param control_network If \code{geneSizeControl=TRUE},
#' then must provide the control network.
#' @inheritParams bootstrap_enrichment_test
#'
#' @returns A list containing three elements:
#' \itemize{
#'   \item \code{hit.cells}: vector containing the summed proportion of
#'   expression in each cell type for the target list.
#'   \item \code{gene_data: } data.table showing the number of time each gene 
#'    appeared in the bootstrap sample.
#'   \item \code{bootstrap_data}: matrix in which each row represents the
#'   summed proportion of expression in each cell type for one of the
#'   random lists
#'   \item \code{controlledCT}: the controlled cell type (if applicable)
#' }
#'
#' @keywords internal
#' @importFrom data.table data.table rbindlist
#' @importFrom parallel mclapply
get_summed_proportions <- function(hits,
                                   sct_data,
                                   annotLevel,
                                   reps,
                                   no_cores = 1,
                                   geneSizeControl,
                                   controlledCT = NULL,
                                   control_network = NULL,
                                   store_gene_data = TRUE,
                                   verbose = TRUE) {
    
    controlledCT <- fix_celltype_names(celltypes = controlledCT)
    combinedGenes <- rownames(sct_data[[annotLevel]]$mean_exp)
    hits <- hits[hits %in% combinedGenes]
    #### cell_list_dist ####
    # cell_list_dist gets the summed proportion of 'hits' across all
    # cell types at annotLevel
    hit.cells <- cell_list_dist(
        hits = hits,
        sct_data = sct_data,
        annotLevel = annotLevel
    )
    #### Check control_network provided if geneSizeControl=TRUE ####
    err_msg <- paste0(
        "ERROR: if geneSizeControl==TRUE then",
        " get_summed_proportions must be passed",
        " control_network as an argument"
    )
    if (isTRUE(geneSizeControl)) {
        if (is.null(control_network)) {
            stop(err_msg)
        }
    }
    # If controlling for expression of another cell type, then the samples of
    # bootstrapped genes must match the deciles of expression specificity
    if (!is.null(controlledCT)) {
        controlled_bootstrap_set <-
            generate_controlled_bootstrap_geneset(
                hits = hits,
                sct_data = sct_data,
                annotLevel = annotLevel,
                reps = reps,
                controlledCT = controlledCT
            )
    }
    #### Parallelise bootstrapping ####
    bootstrap_list <- get_summed_proportions_iterate(
        reps=reps,
        geneSizeControl=geneSizeControl,
        control_network=control_network,
        controlledCT=controlledCT,
        controlled_bootstrap_set=controlled_bootstrap_set,
        combinedGenes=combinedGenes,
        hits=hits,
        sct_data=sct_data,
        annotLevel=annotLevel,
        no_cores=no_cores)
    #### Get gene scores #####
    if(isTRUE(store_gene_data)){
      gene_data <- compute_gene_scores(sct_data = sct_data, 
                                       annotLevel = annotLevel, 
                                       bootstrap_list = bootstrap_list, 
                                       hits = hits, 
                                       combinedGenes = combinedGenes,
                                       verbose = verbose)
    } else {
      gene_data <- NULL
    } 
    #### Get celltypes scores ####
    bootstrap_data <- as.matrix(
        lapply(bootstrap_list,function(x){x$celltypes}) |>
            data.table::rbindlist(),
        rownames = paste0("rep",seq_len(reps)))
    return(list(
        hit.cells = hit.cells,
        gene_data = gene_data,
        bootstrap_data = bootstrap_data,
        controlledCT = controlledCT
    ))
}
