check_controlled_args <- function(bg,
    sct_data,
    annotLevel,
    disease_genes,
    hitGenes,
    functional_genes,
    funcGenes,
    combinedGenes) {
    err_msg <- paste0(
        "ERROR: no bg are present in the single cell",
        " dataset. Perhaps it is from the wrong species?"
    )
    if (sum(bg %in% rownames(sct_data[[annotLevel]]$mean_exp),
        na.rm = TRUE
    ) == 0) {
        stop(err_msg)
    }
    err_msg2 <- paste0(
        "ERROR: no disease_genes are present in the single cell",
        " dataset. Perhaps it is from the wrong species?"
    )
    if (sum(disease_genes %in% combinedGenes,
        na.rm = TRUE
    ) == 0) {
        stop(err_msg2)
    }
    err_msg3 <- paste0(
        "ERROR: insufficient disease_genes. Must provide at",
        " least five that are present in the background",
        " gene set & single cell dataset"
    )
    if (sum(hitGenes %in% combinedGenes, na.rm = TRUE) < 5) {
        stop(err_msg3)
    }
    err_msg4 <- paste0(
        "ERROR: no functional_genes are present in the",
        " single cell dataset. Perhaps it is from the",
        " wrong species?"
    )
    if (sum(functional_genes %in% combinedGenes, na.rm = TRUE) == 0) {
        stop(err_msg4)
    }
    err_msg5 <- paste0(
        "ERROR: insufficient functional_genes Must provide",
        " at least five that are present in the background",
        " gene set & single cell dataset"
    )
    if (sum(funcGenes %in% combinedGenes, na.rm = TRUE) < 5) {
        stop(err_msg5)
    }
}
