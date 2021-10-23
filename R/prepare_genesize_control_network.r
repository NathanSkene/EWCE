#' Prepare genesize control network
#'
#' \code{prepare_genesize_control_network} takes a gene list and finds
#' semi-randomly selected gene lists which are matched for gene length and
#' GC content.
#'
#' @param numBOOT Number of gene lists to sample.
#' @inheritParams bootstrap_enrichment_test
#' @inheritParams orthogene::create_background
#'
#' @return A list containing three data frames:
#' \itemize{
#'   \item \code{hitGenes}: Array of HGNC symbols containing the hit genes.
#'   May be slightly reduced if gene length / GC content could not be found
#'   for all genes.
#'   \item \code{list_network}: The control gene lists as a data frame of HGNC
#'   symbols
#' }
#'
#' @keywords internal
#' @importFrom ewceData ensembl_transcript_lengths_GCcontent
#' @importFrom  stats aggregate quantile
#' @importFrom orthogene create_background
prepare_genesize_control_network <- function(hits,
                                             bg = NULL,
                                             numBOOT = 10000,
                                             sctSpecies = NULL,
                                             genelistSpecies = NULL,
                                             verbose = TRUE) {
    target <- Gene.Symbol <- NULL

    ## Currently, can only convert to human genes
    ## since the gene metadata is only for human genes.
    output_species <- "human"
    #### Check species ####
    species <- check_species(
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies,
        verbose = verbose
    )
    genelistSpecies <- species$genelistSpecies
    sctSpecies <- species$sctSpecies
    #### Check background ####
    bg <- orthogene::create_background(
        species1 = sctSpecies,
        species2 = genelistSpecies,
        output_species = output_species,
        bg = bg,
        verbose = verbose
    )
    #### Prepare to query all_hgnc_wtEnsembl ####
    combined_human_genes <- unique(c(hits, bg))
    #### First get all Ensembl gene IDs from the human genes ####
    ## Old method (can only get human genes)
    # all_hgnc_wtEnsembl <- ewceData::all_hgnc_wtEnsembl()
    ## New method (can get any species' genes)
    # Must use gprofiler, as it's the only method with Ensembl IDs
    all_hgnc_wtEnsembl <- orthogene::all_genes(
        species = output_species,
        method = "gprofiler",
        ensure_filter_nas = FALSE,
        verbose = verbose
    ) %>%
        dplyr::rename(
            ensembl_gene_id = target,
            gene_symbol = Gene.Symbol
        )
    messager(
        formatC(length(unique(all_hgnc_wtEnsembl$ensembl_gene_id)),
            big.mark = ","
        ),
        output_species, "Ensembl IDs", "and",
        formatC(length(unique(all_hgnc_wtEnsembl$gene_symbol)),
            big.mark = ","
        ),
        output_species, "Gene Symbols",
        "imported.",
        v = verbose
    )
    hum_ens <- all_hgnc_wtEnsembl[
        all_hgnc_wtEnsembl$gene_symbol %in% combined_human_genes,
    ]
    #### Check the gene lists are proper gene symbols ####
    if (sum(hits %in% hum_ens$gene_symbol) == 0) {
        err_msg <- paste0(
            "ERROR: No hits recognised as human HGNC symbols.",
            " Perhaps the gene list was wrongly provided as MGI",
            " symbols? prepare_genesize_control_network only",
            " accepts HGNC symbols."
        )
        stop(err_msg)
    }
    if (sum(bg %in% hum_ens$gene_symbol) == 0) {
        err_msg2 <- paste0(
            "ERROR: No bg genes recognised as human HGNC symbols.",
            " Perhaps the gene list was wrongly provided as MGI",
            " symbols? prepare_genesize_control_network only",
            " accepts HGNC symbols."
        )
        stop(err_msg2)
    }
    ### GET THE TRANSCRIPT LENGTHS AND GC CONTENT FROM ewceData ####
    ensembl_transcript_lengths_GCcontent <-
        ewceData::ensembl_transcript_lengths_GCcontent()
    all_lengths <-
        ensembl_transcript_lengths_GCcontent[
            ensembl_transcript_lengths_GCcontent$ensembl_gene_id %in%
                hum_ens$ensembl_gene_id,
        ]
    all_lengths <- all_lengths[!is.na(all_lengths$transcript_length), ]
    all_lens <- merge(
        x = all_lengths,
        y = hum_ens,
        by = "ensembl_gene_id"
    )
    #### TAKE THE MEAN TRANSCRIPT LENGTH & GC-CONTENT FOR EACH GENE ####
    transcript_lengths <- stats::aggregate(
        x = all_lens$transcript_length,
        by = list(all_lens$gene_symbol),
        FUN = mean
    )
    percentage_gene_gc_content <- stats::aggregate(
        x = all_lens$percentage_gene_gc_content,
        by = list(all_lens$gene_symbol),
        FUN = mean
    )
    data_byGene <- cbind(
        transcript_lengths,
        percentage_gene_gc_content[, 2]
    )
    colnames(data_byGene) <- c(
        "HGNC.symbol", "transcript_lengths",
        "percentage_gene_gc_content"
    )
    data_byGene <- data_byGene[data_byGene$HGNC.symbol != "", ]
    #### Drop genes not in 1:1 ortholog-derived bg ###
    # Previously used less flexible ewceData::mouse_to_human_homologs()
    data_byGene2 <- data_byGene[
        data_byGene$HGNC.symbol %in% combined_human_genes,
    ]
    #### Assign each gene to gene length/GC content quadrant ####
    data_byGene2 <- create_quadrants(data_byGene2 = data_byGene2)
    ## FOR EACH 'DISEASE LIST' GENERATE SET
    ## OF 10000 QUADRANT MATCHED GENE LISTS
    ## -Get new set of sctSpecies hitGenes,
    ## containing only those within data_byGene2
    hitGenes_NEW <- data_byGene2[
        data_byGene2$HGNC.symbol %in% hits,
        "HGNC.symbol"
    ]
    list_genes1d <- hits[hits %in% data_byGene$HGNC.symbol]
    ## Get matrix with 10000 randomly sampled genes from the same quadrant
    ## as the gene list.
    list_network <- create_list_network(
        data_byGene2 = data_byGene2,
        hitGenes_NEW = hitGenes_NEW,
        numBOOT = numBOOT
    )
    messager("Controlled bootstrapping network generated.", v = verbose)
    return(list(
        hitGenes = hitGenes_NEW,
        list_network = list_network
    ))
}
