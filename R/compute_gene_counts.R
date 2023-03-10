#' Compute gene counts
#' 
#' Counts the number of times each gene appeared in 
#' the randomly sampled gene lists.
#' @param bootstrap_list The output of \code{get_summed_proportions_iterate}.
#' @inheritParams get_summed_proportions
#' @inheritParams bootstrap_enrichment_test
#' @returns \link[data.table]{data.table}
#' 
#' @keywords internal
#' @importFrom data.table data.table as.data.table setorderv :=
compute_gene_counts<- function(bootstrap_list,
                                # hits,
                                verbose=TRUE){
    
    count <- gene <- proportion_reps <- is_hit_gene <- reps <- NULL;
    
    messager("Computing gene counts.",v=verbose)
    gene_counts <- lapply(bootstrap_list,function(x){x$genes}) |>
        unlist() |> table()
    gene_agg <- data.table::as.data.table(gene_counts) |> 
        `colnames<-`(c("gene","count")) 
    gene_agg[,reps:=length(bootstrap_list)]
    #### Add genes that didn't appear in any randomly sampled list ####
    # extra_genes <- unname(combinedGenes[!combinedGenes %in% gene_agg$gene])
    # if(length(extra_genes)>0){
    #     gene_agg <- rbind(
    #         gene_agg,
    #         data.table::data.table(
    #             gene=extra_genes,
    #             count=0,
    #             reps=length(bootstrap_list))
    #     ) 
    # } 
    
    gene_agg[,proportion_reps:=count/reps,] 
    # gene_agg[,is_hit_gene:=gene %in% hits] 
    # data.table::setorderv(gene_agg,
    #                       cols = c("is_hit_gene","proportion_reps"),
    #                       order = c(-1,-1))
    return(gene_agg)
}
