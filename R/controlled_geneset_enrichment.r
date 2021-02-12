#' Celltype controlled geneset enrichment
#'
#' \code{controlled_geneset_enrichment} tests whether a functional geneset is 
#' still enriched in a disease gene set after controlling for the
#' disease geneset's enrichment in a particular cell type (the 'controlledCT')
#'
#' @param sct_data List generated using \code{\link{generate.celltype.data}}
#' @param disease_genes Array of gene symbols containing the disease gene list. 
#' Does not have to be disease genes. Must be from same species as the single 
#' cell transcriptome dataset.
#' @param functional_genes Array of gene symbols containing the functional gene 
#' list. The enrichment of this geneset within the disease_genes is tested. 
#' Must be from same species as the single cell transcriptome dataset.
#' @param bg_genes Array of gene symbols containing the background gene list.
#' @param reps Number of random gene lists to generate (default=100 but should 
#' be over 10000 for publication quality results)
#' @param annotLevel an integer indicating which level of the annotation to 
#' analyse. Default = 1.
#' @param controlledCT (optional) If not NULL, and instead is the name of a 
#' cell type, then the bootstrapping controls for expression within that 
#' cell type
#' @return A list containing three data frames:
#' \itemize{
#'   \item \code{p_controlled} The probability that functional_genes are 
#'   enriched in disease_genes while controlling for the level of specificity 
#'   in controlledCT
#'   \item \code{z_controlled} The z-score that functional_genes are enriched 
#'   in disease_genes while controlling for the level of specificity in 
#'   controlledCT
#'   \item \code{p_uncontrolled} The probability that functional_genes are 
#'   enriched in disease_genes WITHOUT controlling for the level of 
#'   specificity in controlledCT
#'   \item \code{z_uncontrolled} The z-score that functional_genes are enriched 
#'   in disease_genes WITHOUT controlling for the level of specificity in 
#'   controlledCT
#'   \item \code{reps=reps}
#'   \item \code{controlledCT}
#'   \item \code{actualOverlap=actual} The number of genes that overlap between 
#'   functional and disease gene sets
#' }
#' @examples
#' # See the vignette
#' @export
controlled_geneset_enrichment <- function(disease_genes, functional_genes, 
                                            bg_genes, sct_data, annotLevel, 
                                                reps, controlledCT) {
    combinedGenes <- rownames(sct_data[[annotLevel]]$mean_exp)
    combinedGenes <- combinedGenes[combinedGenes %in% bg_genes]
    hitGenes <- disease_genes[disease_genes %in% combinedGenes]
    funcGenes <- functional_genes[functional_genes %in% combinedGenes]
    err_msg <- paste0("ERROR: no bg_genes are present in the single cell",
                        " dataset. Perhaps it is from the wrong species?")
    if (sum(bg_genes %in% rownames(sct_data[[annotLevel]]$mean_exp)) == 0) {
        stop(err_msg)
    }
    err_msg2 <- paste0("ERROR: no disease_genes are present in the single cell",
                        " dataset. Perhaps it is from the wrong species?")
    if (sum(disease_genes %in% combinedGenes) == 0) {
        stop(err_msg2)
    }
    err_msg3 <- paste0("ERROR: insufficient disease_genes. Must provide at",
                        " least five that are present in the background",
                        " gene set & single cell dataset")
    if (sum(hitGenes %in% combinedGenes) < 5) {
        stop(err_msg3)
    }
    err_msg4 <- paste0("ERROR: no functional_genes are present in the",
                        " single cell dataset. Perhaps it is from the",
                        " wrong species?")
    if (sum(functional_genes %in% combinedGenes) == 0) {
        stop(err_msg4)
    }
    err_msg5 <- paste0("ERROR: insufficient functional_genes Must provide",
                        " at least five that are present in the background",
                        " gene set & single cell dataset")
    if (sum(funcGenes %in% combinedGenes) < 5) {
        stop(err_msg5)
    }

    # Drop genes from sctData which are not in the background set
    sct_data_bgReduced <- sct_data
    for (lvl in seq_len(length(sct_data_bgReduced))) {
        sct_data_bgReduced[[lvl]]$mean_exp <- 
            sct_data_bgReduced[[lvl]]$mean_exp[combinedGenes, ]
        sct_data_bgReduced[[lvl]]$specificity <- 
            sct_data_bgReduced[[lvl]]$specificity[combinedGenes, ]
    }

    controlled_bootstrap_set <- 
        generate_controlled_bootstrap_geneset(hitGenes = hitGenes, 
                                                sct_data = sct_data_bgReduced, 
                                                annotLevel = annotLevel, 
                                                reps = reps, 
                                                controlledCT = controlledCT)
    uncontrolled_bootstrap_set <- 
        replicate(reps, sample(combinedGenes, length(hitGenes)))

    # Calculate uncontrolled enrichment
    numInList <- function(x) {
        return(sum(funcGenes %in% x))
    }
    actual <- numInList(hitGenes)
    bootstrap_uncontrolled <- apply(uncontrolled_bootstrap_set, 2, numInList)
    p_uncontrolled <- sum(bootstrap_uncontrolled > actual) / reps
    z_uncontrolled <- 
        (actual - mean(bootstrap_uncontrolled)) / sd(bootstrap_uncontrolled)

    # Calculate controlled enrichment
    bootstrap_controlled <- apply(controlled_bootstrap_set, 2, numInList)
    p_controlled <- sum(bootstrap_controlled > actual) / reps
    z_controlled <- 
        (actual - mean(bootstrap_controlled)) / sd(bootstrap_controlled)

    # Return results
    return(list(p_controlled = p_controlled, z_controlled = z_controlled, 
                p_uncontrolled = p_uncontrolled, 
                z_uncontrolled = z_uncontrolled, reps = reps, 
                controlledCT = controlledCT, actualOverlap = actual))
}
