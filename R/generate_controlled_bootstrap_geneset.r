#' generate_controlled_bootstrap_geneset
#'
#' \code{generate_controlled_bootstrap_geneset} Used to generated celltype controlled bootstraped
#'
#' @param hitGenes Array of gene names. The target gene set.
#' @param sct_data The cell type data list (with specificity and mean_exp)
#' @param annotLevel The level of annotation in sct_data to analyse
#' @param reps The number of gene lists to sample
#' @param controlledCT Name of a celltype (from colnames of sct_data[[x]]$specificity.
#' @examples
#' # See vignette
#' @export
#' @import stats
generate_controlled_bootstrap_geneset <- function(hitGenes,sct_data,annotLevel,reps,controlledCT=NULL){
    if(is.null(controlledCT)){
        stop("ERROR: controlledCT cannot be NULL in generate_controlled_bootstrap_geneset")
    }
    if(annotLevel>length(sct_data)){
        stop("ERROR: annotLevel cannot be greater than the number of annotation levels in sct_data")
    }

    combinedGenes = rownames(sct_data[[annotLevel]]$mean_exp)
    hitGenes = hitGenes[hitGenes %in% combinedGenes]
    if(length(hitGenes)==0){
        stop("ERROR: length(hitGenes)==0. Perhaps your gene list is from the wrong species? It should be converted to orthologs of the same species as the single cell dataset")
    }
    hit.cells = cell.list.dist(hitGenes,sct_data,annotLevel) # cell.list.dist gets the summed proportion of 'hitGenes' across all cell types at annotLevel

    deciles = unique(quantile(sct_data[[annotLevel]]$specificity[,controlledCT],probs=seq(from=0,to=1,by=0.001)))
    #deciles = unique(quantile(sct_data[[annotLevel]]$specificity[,controlledCT],probs=seq(from=0,to=1,by=0.1)))
    minCount = 0
    for(decile in 1:(length(deciles)-1)){
        decile_min = deciles[decile]
        decile_max = deciles[decile+1]

        # How many of the hitGenes are in this decile?
        hitGenes_specificity = sct_data[[annotLevel]]$specificity[hitGenes,controlledCT]
        num_hitGenes_in_decile = sum(hitGenes_specificity>=decile_min & hitGenes_specificity<decile_max)
        if(num_hitGenes_in_decile>0){

            # Find all the genes in this decile
            allGenes_specificity = sct_data[[annotLevel]]$specificity[combinedGenes,controlledCT]
            genes_in_decile = names(allGenes_specificity)[allGenes_specificity>=decile_min & allGenes_specificity<decile_max]
            decile_boot = replicate(reps,sample(genes_in_decile,num_hitGenes_in_decile))
            minCount = minCount+1
            if(minCount==1){
                controlled_bootstrap_set = decile_boot
            }else{
                controlled_bootstrap_set = rbind(controlled_bootstrap_set,decile_boot)
            }

        }
    }
    return(controlled_bootstrap_set)
}
