#' Celltype controlled geneset enrichment
#'
#' \code{controlled_geneset_enrichment} tests whether a functional gene set is
#' still enriched in a disease gene set after controlling for the
#' disease gene set's enrichment in a particular cell type (the 'controlledCT')
#'
#' @param disease_genes_species Species of the
#' \code{disease_genes} gene set.
#' @param functional_genes_species Species of the
#' \code{functional_genes} gene set.
#' @param disease_genes Array of gene symbols containing the disease gene list.
#' Does not have to be disease genes. Must be from same species as the single
#' cell transcriptome dataset.
#' @param functional_genes Array of gene symbols containing the functional gene
#' list. The enrichment of this gene set within the disease_genes is tested.
#' Must be from same species as the single cell transcriptome dataset.
#' @inheritParams bootstrap_enrichment_test
#' @inheritParams orthogene::create_background
#'
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
#' # See the vignette for more detailed explanations
#' # Gene set enrichment analysis controlling for cell type expression
#' # set seed for bootstrap reproducibility
#' set.seed(12345678)
#' ## load merged dataset from vignette
#' ctd <- ewceData::ctd()
#' schiz_genes <- ewceData::schiz_genes()
#' hpsd_genes <- ewceData::hpsd_genes()
#' # Use 3 bootstrap lists for speed, for publishable analysis use >10000
#' reps <- 3
#'
#' res_hpsd_schiz <- EWCE::controlled_geneset_enrichment(
#'     disease_genes = schiz_genes,
#'     functional_genes = hpsd_genes,
#'     sct_data = ctd,
#'     annotLevel = 1,
#'     reps = reps,
#'     controlledCT = "pyramidal CA1"
#' )
#' @export
#' @importFrom stats sd
controlled_geneset_enrichment <- function(disease_genes,
    functional_genes,
    bg = NULL,
    sct_data,
    sctSpecies = NULL,
    output_species = "human",
    disease_genes_species = NULL,
    functional_genes_species = NULL,
    annotLevel,
    reps = 100,
    controlledCT,
    use_intersect = FALSE,
    verbose = TRUE) {
    #### Check species1 ###
    species <- check_species(
        genelistSpecies = disease_genes_species,
        sctSpecies = sctSpecies,
        verbose = verbose
    )
    disease_genes_species <- species$genelistSpecies
    sctSpecies <- species$sctSpecies
    #### Check species2 ###
    species <- check_species(
        genelistSpecies = functional_genes_species,
        sctSpecies = sctSpecies,
        verbose = verbose
    )
    functional_genes_species <- species$genelistSpecies
    sctSpecies <- species$sctSpecies
    #### Create multi-list bg ####
    ### Also converts gene lists into their output_species orthologs ####
    bg_out <- create_background_multilist(
        gene_list1 = disease_genes,
        gene_list2 = functional_genes,
        gene_list1_species = disease_genes_species,
        gene_list2_species = functional_genes_species,
        output_species = output_species,
        bg = bg,
        use_intersect = use_intersect,
        verbose = verbose
    )
    #### Use background generated in create_background_multilist ####
    bg <- bg_out$bg
    #### Replace gene lists with their output_gene orthologs ####
    disease_genes <- unname(bg_out$gene_list1)
    functional_genes <- unname(bg_out$gene_list2)
    #### Standardise CTD ####
    sct_data <- standardise_ctd(
        ctd = sct_data,
        dataset = "sct_data",
        input_species = sctSpecies,
        output_species = output_species,
        verbose = verbose
    )
    sctSpecies <- output_species
    #### Created combinedGenes list ####
    combinedGenes <- unname(rownames(sct_data[[annotLevel]]$mean_exp))
    combinedGenes <- combinedGenes[combinedGenes %in% bg]
    hitGenes <- disease_genes[disease_genes %in% combinedGenes]
    funcGenes <- functional_genes[functional_genes %in% combinedGenes]
    #### Check args ####
    check_controlled_args(
        bg = bg,
        sct_data = sct_data,
        annotLevel = annotLevel,
        disease_genes = disease_genes,
        hitGenes = hitGenes,
        functional_genes = functional_genes,
        funcGenes = funcGenes,
        combinedGenes = combinedGenes
    )
    #### Drop genes from sctData which are not in the background set ####
    sct_data_bgReduced <- sct_data
    for (lvl in seq_len(length(sct_data_bgReduced))) {
        sct_data_bgReduced[[lvl]]$mean_exp <-
            sct_data_bgReduced[[lvl]]$mean_exp[combinedGenes, ]
        sct_data_bgReduced[[lvl]]$specificity <-
            sct_data_bgReduced[[lvl]]$specificity[combinedGenes, ]
    }
    #### Generate controlled bootstrap gene sets ####
    controlled_bootstrap_set <-
        generate_controlled_bootstrap_geneset(
            hitGenes = hitGenes,
            sct_data = sct_data_bgReduced,
            combinedGenes = combinedGenes,
            annotLevel = annotLevel,
            reps = reps,
            controlledCT = controlledCT,
            verbose = verbose
        )
    uncontrolled_bootstrap_set <- replicate(
        n = reps,
        expr = sample(
            combinedGenes,
            length(hitGenes)
        )
    )
    #### Calculate uncontrolled enrichment ####
    numInList <- function(x) {
        return(sum(funcGenes %in% x))
    }
    actual <- numInList(hitGenes)
    bootstrap_uncontrolled <- apply(uncontrolled_bootstrap_set, 2, numInList)
    p_uncontrolled <- sum(bootstrap_uncontrolled > actual) / reps
    z_uncontrolled <-
        (actual - mean(bootstrap_uncontrolled)) /
            stats::sd(bootstrap_uncontrolled)
    #### Calculate controlled enrichment ####
    bootstrap_controlled <- apply(controlled_bootstrap_set, 2, numInList)
    p_controlled <- sum(bootstrap_controlled > actual) / reps
    z_controlled <-
        (actual - mean(bootstrap_controlled)) /
            stats::sd(bootstrap_controlled)
    #### Return results ####
    return(list(
        p_controlled = p_controlled,
        z_controlled = z_controlled,
        p_uncontrolled = p_uncontrolled,
        z_uncontrolled = z_uncontrolled,
        reps = reps,
        controlledCT = controlledCT,
        actualOverlap = actual
    ))
}
