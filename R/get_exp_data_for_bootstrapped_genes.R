get_exp_data_for_bootstrapped_genes <- function(results,
    signif_res,
    sct_data,
    mouse.hits,
    combinedGenes,
    annotLevel,
    nReps = NA) {
    exp_mats <- list()
    for (cc in signif_res) {
        exp_mats[[cc]] <- matrix(0, nrow = nReps, ncol = length(mouse.hits))
        rownames(exp_mats[[cc]]) <- sprintf("Rep%s", seq_len(nReps))
    }
    for (s in seq_len(nReps)) {
        bootstrap_set <- sample(combinedGenes, length(mouse.hits))
        ValidGenes <- rownames(sct_data[[annotLevel]]$specificity)[
            rownames(sct_data[[annotLevel]]$specificity) %in% bootstrap_set
        ]

        expD <- sct_data[[annotLevel]]$specificity[ValidGenes, ]

        for (cc in signif_res) {
            exp_mats[[cc]][s, ] <- sort(expD[, cc])
        }
    }
    return(exp_mats)
}
