#' check_ewce_genelist_inputs
#'
#' \code{check_ewce_genelist_inputs} Is used to check that hits and bg gene
#' lists passed to EWCE are setup correctly. Checks they are the
#' appropriate length.
#' Checks all hits are in bg. Checks the species match and if not
#' reduces to 1:1 orthologs.
#' @param min_genes Minimum number of genes in a gene list to test.
#' @inheritParams bootstrap_enrichment_test
#' @inheritParams orthogene::convert_orthologs
#'
#' @return A list containing
#' \itemize{
#'   \item \code{hits}: Array of MGI/HGNC gene symbols containing the target
#'   gene list.
#'   \item \code{bg}: Array of MGI/HGNC gene symbols containing the background
#'   gene list.
#' }
#'
#' @examples
#' ctd <- ewceData::ctd()
#' example_genelist <- ewceData::example_genelist()
#'
#' # Called from "bootstrap_enrichment_test()" and "generate_bootstrap_plots()"
#' checkedLists <- EWCE::check_ewce_genelist_inputs(
#'     sct_data = ctd,
#'     hits = example_genelist,
#'     sctSpecies = "mouse",
#'     genelistSpecies = "human"
#' )
#' @export
#' @importFrom orthogene create_background map_genes map_orthologs
check_ewce_genelist_inputs <- function(sct_data,
                                       hits,
                                       bg = NULL,
                                       genelistSpecies = NULL,
                                       sctSpecies = NULL,
                                       sctSpecies_origin = sctSpecies,
                                       output_species = "human",
                                       method = "homologene",
                                       geneSizeControl = FALSE,
                                       standardise_sct_data = TRUE,
                                       standardise_hits = FALSE,
                                       min_genes = 4,
                                       verbose = TRUE) {
    messager("Checking gene list inputs.", v = verbose)
    #### Check species ####
    species <- check_species(
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies,
        sctSpecies_origin = sctSpecies_origin,
        verbose = verbose
    )
    genelistSpecies <- species$genelistSpecies
    sctSpecies <- species$sctSpecies
    sctSpecies_origin <- species$sctSpecies_origin
    # geneSizeControl assumes the genesets are from human genetics...
    # so genelistSpecies must equal "human"
    if (isTRUE(geneSizeControl) &
        (genelistSpecies != "human")) {
        err_msg6 <- paste0(
            "geneSizeControl assumes the genesets are from",
            " human genetics... so genelistSpecies must be set to",
            " 'human'"
        )
        stop(err_msg6)
    }
    #### Create background if none provided ####
    # Keep internal bc it has a check_species beforehand
    if(is.null(bg)){
        bg <- orthogene::create_background(
            species1 = sctSpecies_origin,
            species2 = genelistSpecies,
            output_species = output_species,
            bg = bg,
            method = method,
            verbose = verbose
        )
    }
    #### Convert CTD ####
    #### Standardise sct_data ####
    if ((sctSpecies != output_species) &&
        isTRUE(standardise_sct_data)) {
        messager("Standardising sct_data.", v = verbose)
        sct_data <- standardise_ctd(
            ctd = sct_data,
            input_species = sctSpecies_origin,
            output_species = output_species,
            dataset = "sct_data",
            method = method,
            verbose = FALSE
        )
    }
    sct_genes <- unname(rownames(sct_data[[1]]$mean_exp))
    ##### Standardise hits #### 
    hits <- unique(as.character(hits))
    ## Within-species
    if (genelistSpecies == output_species) {
        if(isTRUE(standardise_hits)){
            messager("Converting gene list input to standardised",
                     output_species, "genes.",
                     v = verbose
            )
            hits <- orthogene::map_genes(
                genes = hits,
                species = genelistSpecies,
                drop_na = TRUE,
                verbose = FALSE
            )$name
        } 
    } else {
    ## Across-species
        messager("Converting gene list input to standardised",
                 output_species, "genes.",
                 v = verbose)
        hits <- orthogene::map_orthologs(
            genes = hits,
            input_species = genelistSpecies,
            output_species = output_species,
            method = method,
            standardise_genes = standardise_hits,
            verbose = FALSE
        )$ortholog_gene
    }
    #### Check that all 'hits' are in 'bg' ####
    hits <- hits[hits %in% bg]
    if (sum(!hits %in% bg, na.rm = TRUE) > 0) {
        stop("All hits must be in bg.")
    }
    #### Check that all 'sct_genes' are in 'bg' ####
    sct_genes <- sct_genes[sct_genes %in% bg]
    if (sum(!sct_genes %in% bg, na.rm = TRUE) > 0) {
        stop("All sct_data genes must be in bg.")
    }
    #### Check that sufficient genes are still present in the target list ####
    if(min_genes<1){
        messager("min_genes must be >0. Setting to 1.",v=verbose)
        min_genes <- 1
    }
    if (length(hits) < min_genes) {
        err_msg5 <- paste0(
            "At least ",min_genes," genes which are present in the",
            " single cell dataset & background gene set are",
            " required to test for enrichment.",
            " Only ",length(hits)," provided."
        )
        stop(err_msg5)
    }
    #### Restrict gene sets to only genes in the SCT dataset  ####
    # if (!geneSizeControl) {
        hits <- hits[hits %in% sct_genes]
        bg <- bg[bg %in% sct_genes]
    # }
    #### Remove all hit genes from bg ####
    bg <- bg[!bg %in% hits]
    #### Return list ####
    return(list(
        hits = hits,
        sct_genes = sct_genes,
        sct_data = sct_data,
        bg = bg,
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies,
        output_species = output_species
    ))
}
