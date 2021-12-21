#' get_exp_data_for_bootstrapped_genes
#' 
#' Support function for 
#' \link[EWCE]{generate_bootstrap_plots_for_transcriptome}.
#' 
#' @param full_results full_results (#fix).
#' @param signif_res signif_res (#fix).
#' @param hits Gene hits.
#' @param combinedGenes Combined list of genes from \code{sct_data}, 
#' \code{hits}, and background \code{bg}. 
#' @inheritParams generate_bootstrap_plots_for_transcriptome
#' @returns exp_mats
#' 
#' @keywords internal
get_exp_data_for_bootstrapped_genes <- function(results,
                                                signif_res,
                                                sct_data,
                                                hits,
                                                combinedGenes,
                                                annotLevel, 
                                                nReps = 100,
                                                as_sparse = TRUE,
                                                verbose = TRUE) {
    messager("Generating exp data for bootstrap genes.",v=verbose)
    #### Extract specificity matrix ####
    spec <- sct_data[[annotLevel]]$specificity
    sct_genes <- rownames(spec)
    #### intialize empty matrices ####
    exp_mats <- list()
    for(cc in signif_res){
        exp_mats[[cc]] <- matrix(0, nrow = nReps, ncol = length(hits))
        rownames(exp_mats[[cc]]) <- sprintf("Rep%s", seq_len(nReps)) 
    }
    #### populate matrices ####
   
    for(s in seq_len(nReps)){
        bootstrap_set <- sample(combinedGenes, length(hits))
        ValidGenes <- sct_genes[sct_genes %in% bootstrap_set] 
        expD <- spec[ValidGenes, ] 
        for (cc in signif_res) {
            exp_mats[[cc]][s, ] <- sort(expD[, cc])
        } 
    }
    #### Convert to sparse matrices ####
    if(as_sparse){
        messager("Converting data for bootstrap tests to sparse matrices.",
                 v=verbose)
        for(cc in signif_res){
            exp_mats[[cc]] <- to_sparse_matrix(exp = exp_mats[[cc]], 
                                               verbose = FALSE)
        }
    } 
    return(exp_mats)
}
