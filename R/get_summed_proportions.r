#' Get summed proportions
#'
#' \code{get_summed_proportions} Given the target geneset, randomly sample 
#' genelists of equal length, obtain the specificity of these and then
#' obtain the mean specificity in each sampled list (and the target list)
#'
#' @param hitGenes Array of gene names. The target gene set.
#' @param sct_data The cell type data list (with specificity and mean_exp)
#' @param annotLevel The level of annotation in sct_data to analyse
#' @param reps The number of gene lists to sample
#' @param geneSizeControl Boolean. Should the sampled gene lists control for 
#' GC content and transcript length?
#' @param controlledCT Name of a celltype (from colnames of 
#' sct_data[[x]]$specificity).
#' @param control_network If geneSizeControl=TRUE, then must provide the 
#' control network
#' @return A list containing three data frames:
#' \itemize{
#'   \item \code{hit.cells}: vector containing the summed proportion of 
#'   expression in each cell type for the target list
#'   \item \code{bootstrap_data}: matrix in which each row represents the 
#'   summed proportion of expression in each cell type for one of the 
#'   random lists
#'   \item \code{controlledCT}: the controlled celltype (if applicable)
#' }
#' @examples
#' # See bootstrap_enrichment_test.r
get_summed_proportions <- function(hitGenes, sct_data, annotLevel, reps, 
                                    geneSizeControl, controlledCT = NULL, 
                                    control_network = NULL) {
    bootstrap_data <- matrix(0, ncol = length(
            colnames(sct_data[[annotLevel]]$specificity)),
            nrow = reps)
    combinedGenes <- rownames(sct_data[[annotLevel]]$mean_exp)
    hitGenes <- hitGenes[hitGenes %in% combinedGenes]
    hit.cells <- cell_list_dist(hitGenes, sct_data, annotLevel)
    # cell_list_dist gets the summed proportion of 'hitGenes' across all 
    # cell types at annotLevel

    # Check control_network provided if geneSizeControl=TRUE
    err_msg <- paste0("ERROR: if geneSizeControl==TRUE then",
                        " get_summed_proportions must be passed",
                        " control_network as an argument")
    if (geneSizeControl == TRUE) {
        if (is.null(control_network)) {
            stop(err_msg)
        }
    }

    # If controlling for expression of another cell type, then the samples of 
    # bootstrapped genes must match the deciles of expression specificity
    if (!is.null(controlledCT)) {
        controlled_bootstrap_set <- 
            generate_controlled_bootstrap_geneset(hitGenes, sct_data, 
                                                    annotLevel, reps, 
                                                    controlledCT)
    }

    for (s in seq_len(reps)) {
        # Get 'bootstrap_set'...a list of genes of equivilent length as hitGenes
        if (geneSizeControl == FALSE) {
            if (is.null(controlledCT)) {
                bootstrap_set <- sample(combinedGenes, length(hitGenes))
            } else {
                bootstrap_set <- controlled_bootstrap_set[, s]
            }
        } else {
            bootstrap_set <- control_network[s, ]
        }
        # 'bootstrap_data' is a matrix of the summed proportions
        bootstrap_data[s, ] <- cell_list_dist(bootstrap_set, sct_data, 
                                                annotLevel)
    }
    colnames(bootstrap_data) <- names(hit.cells)

    return(list(hit.cells = hit.cells, bootstrap_data = bootstrap_data, 
                    controlledCT = controlledCT))
}
