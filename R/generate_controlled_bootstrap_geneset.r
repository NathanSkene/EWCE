#' generate_controlled_bootstrap_geneset
#'
#' Used to generated cell type-controlled bootstrapped gene sets.
#'
#' See \link[EWCE]{controlled_geneset_enrichment} for examples.
#'
#' @param hitGenes Array of gene names. The target gene set.
#' @param sct_data The cell type data list (with specificity and mean_exp)
#' @param annotLevel The level of annotation in sct_data to analyse
#' @param reps The number of gene lists to sample
#' @param controlledCT Name of a cell type (from colnames of
#' \code{sct_data[[x]]$specificity}).
#' @param verbose Print messages.
#'
#' @return Matrix of genes
#'  (such that \code{nrows=length(hitGenes)} and \code{ncols=reps}), where each
#' column is a gene list.
#'
#' @keywords internal
#' @importFrom stats quantile
generate_controlled_bootstrap_geneset <- function(hitGenes,
                                                  sct_data,
                                                  combinedGenes,
                                                  annotLevel,
                                                  reps,
                                                  controlledCT = FALSE,
                                                  verbose = TRUE) {
    messager("Generating controlled bootstrap gene sets.", v = verbose)
    err_msg <- paste0(
        "ERROR: controlledCT cannot be NULL in",
        " generate_controlled_bootstrap_geneset"
    )
    if (is.null(controlledCT)) {
        stop(err_msg)
    }
    err_msg2 <- paste0(
        "ERROR: annotLevel cannot be greater than the number",
        " of annotation levels in sct_data"
    )
    if (annotLevel > length(sct_data)) {
        stop(err_msg2)
    }
    # Check all controlledCT are in single cell data
    err_msg3 <- paste0(
        "ERROR: not all controlledCT are in",
        " colnames(sct_data[[annotLevel]]$specificity)"
    )
    if (sum(!controlledCT %in%
        colnames(sct_data[[annotLevel]]$specificity)) != 0) {
        stop(err_msg3)
    }
    err_msg4 <- paste0(
        "ERROR: length(hitGenes)==0. Perhaps your gene list is",
        " from the wrong species? It should be converted to",
        " orthologs of the same species as the single cell",
        " dataset"
    )
    if (length(hitGenes) == 0) {
        stop(err_msg4)
    }
    hit.cells <- cell_list_dist(
        hitGenes = hitGenes,
        sct_data = sct_data,
        annotLevel = annotLevel
    )
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
        tmp_deciles <- stats::quantile(
            sct_data[[annotLevel]]$specificity[, cCT],
            probs = quantile_probs
        )
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

    # For each gene, find its specificity in each controlled cell type
    eachGeneSP <- matrix(0,
        nrow = dim(sct_data[[annotLevel]]$specificity)[1],
        ncol = length(controlledCT)
    )
    rownames(eachGeneSP) <- rownames(sct_data[[annotLevel]]$specificity)
    colnames(eachGeneSP) <- controlledCT
    for (cCT in controlledCT) {
        for (gg in rownames(eachGeneSP)) {
            geneSpecificity <- sct_data[[annotLevel]]$specificity[gg, cCT]
            whichIDX <- sort(which(ct_deciles[, cCT] < geneSpecificity),
                decreasing = TRUE
            )[1]
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
            controlled_bootstrap_set <- rbind(
                controlled_bootstrap_set,
                decile_boot
            )
        }
    }

    return(controlled_bootstrap_set)
}
