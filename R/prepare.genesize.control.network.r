#' Prepare genesize control network
#'
#' \code{prepare.genesize.control.network} takes a genelist and finds 
#' semi-randomly selected gene lists which are matched for gene length and 
#' GC content
#'
#' @param hits Array of MGI gene symbols containing the target gene list. 
#' Must be HGNC symbols.
#' @param bg Array of MGI gene symbols containing the background gene list 
#' (including hit genes). Must be HGNC symbols.
#' @param numBOOT Number of gene lists to sample
#' @param sctSpecies  Either 'mouse' or 'human' depending on whether MGI or 
#' HGNC symbols are used for the single cell dataset
#' @return A list containing three data frames:
#' \itemize{
#'   \item \code{hitGenes}: Array of HGNC symbols containing the hit genes. 
#'   May be slightly reduced if gene length / GC content could not be found 
#'   for all genes.
#'   \item \code{list_network}: The control gene lists as a data frame of HGNC 
#'   symbols
#' }
#' @examples
#' # Called by bootstrap.enrichment.t
#' @import stats
#' @import utils
#' @import biomaRt
#' @import ewceData
prepare.genesize.control.network <- function(hits, bg, numBOOT = 10000, 
                                                sctSpecies) {
    ### PREPARE TO QUERY BIOMART
    combined_human_genes <- unique(c(hits, bg))

    ### FIRST GET ALL ENSEMBL GENE IDS FOR THE HUMAN GENES
    all_hgnc_wtEnsembl <- ewceData::all_hgnc_wtEnsembl()
    hum_ens <- 
        all_hgnc_wtEnsembl[
            all_hgnc_wtEnsembl$hgnc_symbol %in% 
                combined_human_genes, ]

    # CHECK THE GENELISTS WERE HUMAN HGNC SYMBOLS
    err_msg <- paste0("ERROR: No hits recognised as human HGNC symbols.",
                        " Perhaps the gene list was wrongly provided as MGI",
                        " symbols? prepare.genesize.control.network only",
                        " accepts HGNC symbols.")
    err_msg2 <- paste0("ERROR: No bg genes recognised as human HGNC symbols.",
                        " Perhaps the gene list was wrongly provided as MGI",
                        " symbols? prepare.genesize.control.network only",
                        " accepts HGNC symbols.")
    if (sum(hits %in% hum_ens$hgnc_symbol) == 0) {
        stop(err_msg)
    }
    if (sum(bg %in% hum_ens$hgnc_symbol) == 0) {
        stop(err_msg2)
    }

    ### GET THE TRANSCRIPT LENGTHS AND GC CONTENT FROM BIOMART
    ensembl_transcript_lengths_GCcontent <- 
        ewceData::ensembl_transcript_lengths_GCcontent()
    all_lengths <- 
        ensembl_transcript_lengths_GCcontent[
            ensembl_transcript_lengths_GCcontent$ensembl_gene_id %in% 
                hum_ens$ensembl_gene_id, ]
    all_lengths <- all_lengths[!is.na(all_lengths$transcript_length), ]
    all_lens <- merge(all_lengths, hum_ens, by = "ensembl_gene_id")

    # TAKE THE MEAN TRANSCRIPT LENGTH & GC-CONTENT FOR EACH GENE
    transcript_lengths <- aggregate(all_lens$transcript_length, 
                                        by = list(all_lens$hgnc_symbol), 
                                        FUN = mean)
    percentage_gene_gc_content <- aggregate(all_lens$percentage_gene_gc_content,
                                                by = list(all_lens$hgnc_symbol),
                                                FUN = mean)
    data_byGene <- cbind(transcript_lengths, percentage_gene_gc_content[, 2])
    colnames(data_byGene) <- c("HGNC.symbol", "transcript_lengths",
                                "percentage_gene_gc_content")
    data_byGene <- data_byGene[data_byGene$HGNC.symbol != "", ]

    if (sctSpecies == "mouse") {
        ### DROP ANY HUMAN GENES WITHOUT HOMOLOGOUS MOUSE GENES
        m2h <- unique(ewceData::mouse_to_human_homologs()[, c("HGNC.symbol", 
                                                                "MGI.symbol")])
        data_byGene2 <- 
            data_byGene[data_byGene$HGNC.symbol %in% m2h$HGNC.symbol, ]

        ### MERGE THE LENGTH/GC DATA WITH MOUSE ORTHOLOG DATA
        data_byGene3 <- merge(data_byGene2, m2h, by = "HGNC.symbol")
        data_byGene2 <- data_byGene3
    } else if (sctSpecies == "human") {
        data_byGene2 <- data_byGene3 <- data_byGene
    }

    # GET QUANTILES FOR TRANSCRIPT LENGTH + GC CONTENT
    tl_quants <- quantile(data_byGene2$transcript_length, 
                            probs = seq(0.1, 1, 0.1))
    gc_quants <- quantile(data_byGene2$percentage_gene_gc_content, 
                            probs = seq(0.1, 1, 0.1))

    # ASSIGN EACH GENE TO A QUANTILE QUADRANT
    quadrant <- matrix(0, nrow = dim(data_byGene2)[1], ncol = 2)
    colnames(quadrant) <- c("TL", "GC")
    for (i in seq_len(dim(data_byGene2)[1])) {
        quadrant[i, 1] <- which(data_byGene2[i, 2] < tl_quants)[1]
        quadrant[i, 2] <- which(data_byGene2[i, 3] < gc_quants)[1]
    }
    data_byGene2$uniq_quad <- sprintf("%s_%s", quadrant[, 1], quadrant[, 2])
    uq <- data_byGene2$uniq_quad
    data_byGene2 <- data_byGene2[uq != "2_NA" & uq != "NA_2" & uq != "3_NA", ]

    ###FOR EACH 'DISEASE LIST' GENERATE SET OF 10000 QUADRANT MATCHED GENE LISTS
    # -Get new set of mouse hitGenes, containing only those within data_byGene2
    if (sctSpecies == "mouse") {
        hitGenes_NEW <- 
            data_byGene2[data_byGene2$HGNC.symbol %in% hits, "MGI.symbol"]
    } else if (sctSpecies == "human") {
        hitGenes_NEW <- 
            data_byGene2[data_byGene2$HGNC.symbol %in% hits, "HGNC.symbol"]
    }
    list_genes1d <- hits[hits %in% data_byGene$HGNC.symbol]

    # GET ALL MOUSE GENES IN EACH QUADRANT
    quad_genes <- list()
    for (uq in unique(data_byGene2$uniq_quad)) {
        if (sctSpecies == "mouse") {
            quad_genes[[uq]] <- 
                unique(data_byGene2[data_byGene2$uniq_quad == uq, "MGI.symbol"])
        } else if (sctSpecies == "human") {
            quad_genes[[uq]] <- 
                unique(data_byGene2[data_byGene2$uniq_quad == uq, 
                                        "HGNC.symbol"])
        }
    }

    # GET MATRIX WITH 10000 RANDOMLY SAMPLED GENES FROM THE SAME QUADRANT 
    # AS THE LIST GENE
    list_network <- matrix("", nrow = numBOOT, ncol = length(hitGenes_NEW))
    count <- 0
    for (gene in hitGenes_NEW) {
        count <- count + 1
        if (sctSpecies == "mouse") {
            this_gene_quad <- data_byGene2[data_byGene2$MGI.symbol == gene,
                                            "uniq_quad"][1]
        } else if (sctSpecies == "human") {
            this_gene_quad <- data_byGene2[data_byGene2$HGNC.symbol == gene,
                                            "uniq_quad"][1]
        }
        candidates <- as.vector(unlist(quad_genes[this_gene_quad]))
        list_network[, count] <- sample(candidates, numBOOT, replace = TRUE)
    }
    print("CONTROLLED BOOTSTRAPPING NETWORK GENERATED")
    return(list(hitGenes = hitGenes_NEW, list_network = list_network))
}
