#' generate_controlled_bootstrap_geneset
#'
#' \code{generate_controlled_bootstrap_geneset} Used to generated celltype 
#' controlled bootstraped
#'
#' @param hitGenes Array of gene names. The target gene set.
#' @param sct_data The cell type data list (with specificity and mean_exp)
#' @param annotLevel The level of annotation in sct_data to analyse
#' @param reps The number of gene lists to sample
#' @param controlledCT Name of a celltype (from colnames of 
#' sct_data[[x]]$specificity).
#' @return Matrix of genes (nrows=length(hitGenes),ncols=reps), where each 
#' column is a gene list
#' @examples
#' # See controlled_geneset_enrichment.r
#' @import stats
generate_controlled_bootstrap_geneset <- function(hitGenes, sct_data, 
                                                    annotLevel, reps, 
                                                        controlledCT = NULL) {

    err_msg <- paste0("ERROR: controlledCT cannot be NULL in",
                        " generate_controlled_bootstrap_geneset")
    if (is.null(controlledCT)) {
        stop(err_msg)
    }
    err_msg2 <- paste0("ERROR: annotLevel cannot be greater than the number",
                        " of annotation levels in sct_data")
    if (annotLevel > length(sct_data)) {
        stop(err_msg2)
    }
    # Check all controlledCT are in single cell data
    err_msg3 <- paste0("ERROR: not all controlledCT are in",
                        " colnames(sct_data[[annotLevel]]$specificity)")
    if (sum(!controlledCT %in% 
                colnames(sct_data[[annotLevel]]$specificity)) != 0) {
        stop(err_msg3)
    }

    combinedGenes <- rownames(sct_data[[annotLevel]]$mean_exp)
    hitGenes <- hitGenes[hitGenes %in% combinedGenes]
    err_msg4 <- paste0("ERROR: length(hitGenes)==0. Perhaps your gene list is",
                        " from the wrong species? It should be converted to",
                        " orthologs of the same species as the single cell",
                        " dataset")
    if (length(hitGenes) == 0) {
        stop(err_msg4)
    }
    hit.cells <- cell_list_dist(hitGenes, sct_data, annotLevel) 

    # quantile_probs = seq(from=0,to=1,by=0.001)
    if (length(controlledCT) == 1) {
        byStep <- 0.001
    }
    if (length(controlledCT) == 2) {
        byStep <- 0.01
    }
    if (length(controlledCT) >= 3) {
        byStep <- 0.1
    }
    quantile_probs <- seq(from = 0, to = 1, by = byStep)
    for (cCT in controlledCT) {
        tmp_deciles <- quantile(sct_data[[annotLevel]]$specificity[, cCT], 
                                probs = quantile_probs)
        if (cCT == controlledCT[1]) {
            ct_deciles <- tmp_deciles
        } else {
            ct_deciles <- cbind(ct_deciles, tmp_deciles)
        }
    }
    if (is.null(dim(ct_deciles))) {
        ct_deciles <- t(t(ct_deciles))
    }
    colnames(ct_deciles) <- controlledCT
    ct_deciles <- unique(ct_deciles)
    ct_deciles <- ct_deciles[-dim(ct_deciles)[1], , drop = FALSE]

    # For each gene, find it's specificity in each controlled celltype
    eachGeneSP <- matrix(0, nrow = dim(sct_data[[annotLevel]]$specificity)[1], 
                            ncol = length(controlledCT))
    rownames(eachGeneSP) <- rownames(sct_data[[annotLevel]]$specificity)
    colnames(eachGeneSP) <- controlledCT
    for (cCT in controlledCT) {
        for (gg in rownames(eachGeneSP)) {
            geneSpecificity <- sct_data[[annotLevel]]$specificity[gg, cCT]
            whichIDX <- sort(which(ct_deciles[, cCT] < geneSpecificity), 
                                decreasing = TRUE)[1]
            if (is.na(whichIDX)) {
                whichIDX <- 1
            }
            eachGeneSP[gg, cCT] <- ct_deciles[whichIDX, cCT]
        }
    }
    collapseEntries <- function(x) {
        y <- paste(x, collapse = ",")
        return(y)
    }
    eachGeneBOX <- apply(eachGeneSP, 1, collapseEntries)

    boxes_present <- unique(eachGeneBOX)
    boxes_present_inHits <- table(eachGeneBOX[hitGenes])
    # For each box, sample the number of genes as is present in 
    # hitGenes in that box
    minCount <- 0
    for (i in seq_len(length(boxes_present_inHits))) {
        boxName <- names(boxes_present_inHits[i])
        boxFreqInHits <- boxes_present_inHits[i]
        allGenesInBox <- names(eachGeneBOX[eachGeneBOX == boxName])
        decile_boot <- replicate(reps, sample(allGenesInBox, boxFreqInHits))
        minCount <- minCount + 1
        if (minCount == 1) {
            controlled_bootstrap_set <- decile_boot
        } else {
            controlled_bootstrap_set <- rbind(controlled_bootstrap_set, 
                                                decile_boot)
        }
    }

    return(controlled_bootstrap_set)
}
