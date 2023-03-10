get_summed_proportions_iterate <- function(reps,
                                           geneSizeControl,
                                           control_network,
                                           controlledCT,
                                           controlled_bootstrap_set,
                                           combinedGenes,
                                           hits,
                                           sct_data,
                                           annotLevel,
                                           no_cores){
    
    parallel::mclapply(seq_len(reps), function(s) {
        # Get 'bootstrap_set'...a list of genes of equivalent length as hits
        if (isTRUE(geneSizeControl)) {
            bootstrap_set <- control_network[s, ]
        } else {
            if (is.null(controlledCT)) {
                bootstrap_set <- sample(
                    combinedGenes,
                    length(hits)
                )
            } else {
                bootstrap_set <- controlled_bootstrap_set[, s]
            }
        }
        # 'bootstrap_data' is a matrix of the summed proportions
        bootstrap_res <- cell_list_dist(
            hits = bootstrap_set,
            sct_data = sct_data,
            annotLevel = annotLevel
        )
        
        return(list(celltypes=data.table::data.table(t(bootstrap_res)),
                    genes=bootstrap_set
                    )
        )
    }, mc.cores = no_cores)  
}
