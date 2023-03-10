generate_bootstrap_plots_exp_mats <- function(exp_mats=NULL,
                                              sct_data,
                                              annotLevel,
                                              bootstrap_list=NULL,
                                              reps=NULL,
                                              combinedGenes,
                                              hits,
                                              verbose=TRUE){
    
    if(length(exp_mats)>0) return(exp_mats)
    #### Extract precomputed random gene lists, or create new ones ####
    if(!is.null(bootstrap_list) &&
       !is.null(bootstrap_list[[1]]$genes)){
        messager("Using previously sampled genes.",v=verbose)
        reps <- length(bootstrap_list)    
        boot_genes <- lapply(bootstrap_list, 
                             function(i){unname(i$genes)}) 
    } else {
        messager("Resampling random genes.",v=verbose)
        if(is.null(reps)){
            stp <- paste0("Must supply reps when bootstrap_list is NULL.")
            stop(stp)
        } 
        boot_genes <- lapply(seq_len(reps), 
                             function(i){unname(
                                 sample(combinedGenes,
                                        size = length(hits)
                                 ))})
    }
    celltypes <- colnames(sct_data[[annotLevel]]$specificity)
    ### Make template matrices ####
    exp_mats <- lapply(stats::setNames(celltypes,
                                       celltypes),
                       function(cc){
                           matrix(0,
                                  nrow = reps,
                                  ncol = length(hits)
                           ) |> `rownames<-`(sprintf("Rep%s", seq_len(reps)))
                       }) 
    for (s in seq_len(reps)) {
        bootstrap_set <- sample(x = combinedGenes, 
                                size = length(hits))
        ValidGenes <- rownames(sct_data[[annotLevel]]$specificity)[
            rownames(sct_data[[annotLevel]]$specificity) %in% bootstrap_set
        ]
        expD <- sct_data[[annotLevel]]$specificity[ValidGenes, ]
        for (cc in colnames(expD)) {
            exp_mats[[cc]][s, ] <- sort(expD[, cc])
        }
    }
    # exp_dt <- lapply(seq_len(reps), function(i){
    #     ValidGenes <- rownames(sct_data[[annotLevel]]$specificity)[
    #         rownames(sct_data[[annotLevel]]$specificity) %in% boot_genes[[i]]
    #     ]
    #     expD <- sct_data[[annotLevel]]$specificity[ValidGenes, ]
    #     lapply(colnames(expD), function(cc){
    #         sorted <- sort(expD[, cc], decreasing = TRUE)
    #         data.table::data.table(gene=names(sorted), 
    #                                exp=unname(sorted),
    #                                rep=paste0('rep',i))
    #     }) |> data.table::rbindlist(use.names = TRUE, idcol = "CellType") 
    # }) |> data.table::rbindlist() 
   
    return(exp_mats)
}
