#' check_ewce_genelist_inputs
#'
#' \code{check_ewce_genelist_inputs} Is used to check that hits and bg gene
#' lists passed to EWCE are setup correctly. Checks they are the
#' appropriate length.
#' Checks all hits genes are in bg. Checks the species match and if not
#' reduces to 1:1 orthologs.
#'
#' @param sct_data List generated using \code{\link{generate_celltype_data}}
#' @param hits Array of MGI/HGNC gene symbols containing the target gene list.
#' @param bg Array of MGI/HGNC gene symbols containing the background gene list.
#' @param genelistSpecies Either 'mouse' or 'human' depending on whether MGI
#' or HGNC symbols are used for gene lists
#' @param sctSpecies  Either 'mouse' or 'human' depending on whether MGI
#' or HGNC symbols are used for the single cell dataset
#' @param geneSizeControl Will genelists sampled control for GC content and
#' transcript length? Boolean.
#' @return A list containing
#' \itemize{
#'   \item \code{hits}: Array of MGI/HGNC gene symbols containing the target
#'   gene list.
#'   \item \code{bg}: Array of MGI/HGNC gene symbols containing the background
#'   gene list.
#' }
#'
#'
#' @examples
#' library(ewceData)
#' # Called from "bootstrap_enrichment_test()" and "generate_bootstrap_plots()"
#' ctd <- ctd()
#' example_genelist <- example_genelist()
#' mouse_to_human_homologs <- mouse_to_human_homologs()
#' m2h <- unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])
#' mouse.hits <-
#'      unique(m2h[m2h$HGNC.symbol %in% example_genelist, "MGI.symbol"])
#' #subset mouse.bg for speed but ensure it still contains the hits
#' mouse.bg <- unique(c(m2h$MGI.symbol[1:100],mouse.hits))
#' checkedLists <- check_ewce_genelist_inputs(
#'     sct_data = ctd, hits = mouse.hits,
#'     bg = mouse.bg, genelistSpecies = "mouse", sctSpecies = "mouse"
#' )
#' @export
#' @import utils
#' @import ewceData
#' @import ExperimentHub
#' @importFrom AnnotationHub query

check_ewce_genelist_inputs <- function(sct_data, hits, bg, genelistSpecies,
                                        sctSpecies, geneSizeControl = FALSE) {
    # CHECK THE ARGUMENTS ARE PROPERLY STRUCTURED
    orthologsOnly <- FALSE
    sct_genes <- rownames(sct_data[[1]]$mean_exp)

    ## Check that all 'hits' are in 'bg'
    if (sum(!hits %in% bg) > 0) {
        stop("ERROR: all hits must be in bg")
    }

    err_msg <- paste0("ERROR: genelistSpecies and sctSpecies must be either",
                      " mouse or human. If using data from any other species",
                      " then please convert to mouse/human using only",
                      " 1:1 orthologs.")
    ## Check that gene lists and single cell dataset are either mouse or human
    if (sum(c(genelistSpecies, sctSpecies) %in% c("mouse", "human")) != 2) {
        stop(err_msg)
    }

    all_mgi <- ewceData::all_mgi()
    all_hgnc <- ewceData::all_hgnc()

    err_msg2 <- paste0("ERROR: less than four of the hits genes are MGI",
                        " symbols. They must be provided as correctly",
                        " formatted MGI symbols (or alter genelistSpecies)")
    err_msg2_hgnc <- paste0("ERROR: less than four of the hits genes are HGNC",
                             " symbols. They must be provided as correctly",
                             " formatted HGNC symbols (or alter genelistSpecies)")
    err_msg3 <- paste0("ERROR: more bg genes are HGNC genes than MGI genes.",
                        " Did you provide the correct species?")
    err_msg3_hgnc <- paste0("ERROR: more bg genes are MGI genes than HGNC genes.",
                       " Did you provide the correct species?")

    ## Check that the gene lists are really from the correct species
    if (genelistSpecies == "mouse") {
        if (sum(hits %in% all_mgi) < 4) {
            message(sprintf("Passed: %s", paste(hits, collapse = ", ")))
            stop(err_msg2)
        }
        if (sum(bg %in% all_mgi) < 4) {
            message(sprintf("Passed: %s", paste(hits, collapse = ", ")))
            stop(err_msg2)
        }
        if (sum(bg %in% all_mgi) < sum(bg %in% all_hgnc)) {
            stop(err_msg3)
        }
    }
    if (genelistSpecies == "human") {
        if (sum(hits %in% all_hgnc) < 4) {
            stop(err_msg2_hgnc)
        }
        if (sum(bg %in% all_hgnc) < 4) {
            stop(err_msg2_hgnc)
        }
        if (sum(bg %in% all_hgnc) < sum(bg %in% all_mgi)) {
            stop(err_msg3_hgnc)
        }
    }
    err_msg4 <- paste0("ERROR: fewer single cell dataset genes are recognised",
                        " MGI symbols than HGNC symbols. Did you provide the",
                        " correct species? And set sctSpecies correctly?")
    err_msg4_hgnc <- paste0("ERROR: fewer single cell dataset genes are recognised",
                       " HGNC symbols than MGI symbols. Did you provide the",
                       " correct species? And set sctSpecies correctly?")
    if (sctSpecies == "mouse") {
        if (sum(sct_genes %in% all_mgi) < sum(sct_genes %in% all_hgnc)) {
            stop(err_msg4)
        }
    }
    if (sctSpecies == "human") {
        if (sum(sct_genes %in% all_hgnc) < sum(sct_genes %in% all_mgi)) {
            stop(err_msg4_hgnc)
        }
    }

    mouse_to_human_homologs <- ewceData::mouse_to_human_homologs()

    ## If gene lists and single cell data are from different species...
    # then convert the gene lists to match
    if (geneSizeControl == FALSE) {
        if (sctSpecies == "mouse" & genelistSpecies == "human") {
            hits <- mouse_to_human_homologs$MGI.symbol[
                        mouse_to_human_homologs$HGNC.symbol %in% hits]
            bg <- mouse_to_human_homologs$MGI.symbol[
                        mouse_to_human_homologs$HGNC.symbol %in% bg]
            orthologsOnly <- TRUE
            sct_genes <- sct_genes[
                sct_genes %in% mouse_to_human_homologs$MGI.symbol]
        }
    }
    if (sctSpecies == "human" & genelistSpecies == "mouse") {
        hits <- mouse_to_human_homologs$HGNC.symbol[
                    hits %in% mouse_to_human_homologs$MGI.symbol]
        bg <- mouse_to_human_homologs$HGNC.symbol[
                    bg %in% mouse_to_human_homologs$MGI.symbol]
        orthologsOnly <- TRUE
        sct_genes <- sct_genes[
            sct_genes %in% mouse_to_human_homologs$HGNC.symbol]
    }

    # Restrict genesets to only genes in the SCT dataset
    if (geneSizeControl == FALSE) {
        hits <- hits[hits %in% sct_genes]
        bg <- bg[bg %in% sct_genes]
    }
    err_msg5 <- paste0("ERROR: At least four genes which are present in the",
                        " single cell dataset & background gene set are",
                        " required to test for enrichment")
    # Check that sufficient genes are still present in the target list
    if (length(hits) < 4) {
        stop(err_msg5)
    }

    # Remove all hit genes from bg
    bg <- bg[!bg %in% hits]
    err_msg6 <-paste0("ERROR: geneSizeControl assumes the genesets are from",
                        " human genetics... so genelistSpecies must be set to",
                        " 'human'")
    # geneSizeControl assumes the genesets are from human genetics...
    # so genelistSpecies must equal "human"
    if (geneSizeControl == TRUE & genelistSpecies != "human") {
        stop(err_msg6)
    }
    return(list(hits = hits, bg = bg))
}
