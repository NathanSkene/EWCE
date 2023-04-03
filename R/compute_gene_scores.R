#' Compute gene counts
#' 
#' Aggregate gene-level scores across all bootstrap iterations.
#' \itemize{
#' \item{boot: }{Mean specificity of all genes within a given cell type.}
#' \item{hit: }{Mean specificity of a hit gene within a given cell type.}
#' }
#' @param return_hit_exp Return the expression of each hit gene.
#' @param bootstrap_list The output of \code{get_summed_proportions_iterate}.
#' @inheritParams get_summed_proportions
#' @inheritParams bootstrap_enrichment_test
#' @returns \link[data.table]{data.table}
#' 
#' @keywords internal
#' @importFrom data.table data.table as.data.table setorderv := .N
compute_gene_scores <- function(sct_data, 
                                annotLevel,
                                bootstrap_list=NULL,
                                hits,
                                combinedGenes,
                                reps=NULL,
                                exp_mats=NULL,
                                return_hit_exp=FALSE,
                                verbose=TRUE
                                ){
    rank <- p <- hit <- boot <- .N <- NULL;
    
    messager("Computing gene scores.",v=verbose)
    exp_mats <- generate_bootstrap_plots_exp_mats(exp_mats=exp_mats,
                                                  sct_data=sct_data,
                                                  annotLevel=annotLevel,
                                                  bootstrap_list=bootstrap_list,
                                                  reps=reps,
                                                  combinedGenes=combinedGenes,
                                                  hits=hits,
                                                  verbose=verbose)
    #### Get specificity scores of the hit genes ####
    hit_exp <- sct_data[[annotLevel]]$specificity[hits, , drop=FALSE]
    #### Create plotting data ####
    gene_data <- lapply(stats::setNames(names(exp_mats),
                                        names(exp_mats)),
                  function(cc){ 
        mean_boot_exp <-apply(exp_mats[[cc]], 2, mean)
        hit_exp2 <- sort(hit_exp[, cc])
        hit_exp_names <- rownames(hit_exp)[order(hit_exp[, cc])]
        data.table::data.table(
            gene = hit_exp_names,
            boot = mean_boot_exp*100,
            hit = hit_exp2*100,
            rank = order(hit_exp2, decreasing = TRUE)
        )
    }) |> data.table::rbindlist(use.names = TRUE, idcol = "CellType") 
    data.table::setorderv(gene_data,
                          cols = c("CellType","hit"),
                          order = c(1,-1))
    gene_data[,rank:=factor(rank,
                            levels = rev(unique(sort(rank))),
                            ordered = TRUE)]
    #### Compute p-values #####
    gene_data[,p:=sum(hit<boot)/.N, by=c("CellType","rank")] 
    gene_data[,q:=stats::p.adjust(p,method = "bonf"), by=c("CellType","rank")]
    ##### Compute gene counts ####
    if(!is.null(bootstrap_list)){
        gene_counts <- compute_gene_counts(bootstrap_list = bootstrap_list) 
        gene_data <- data.table::merge.data.table(gene_data,
                                                  gene_counts,
                                                  all.x = TRUE,
                                                  by = "gene") 
    } 
    if(isTRUE(return_hit_exp)){
        return(list(gene_data=gene_data,
                    hit_exp=hit_exp))    
    } else{
        return(gene_data)
    }
    
}
