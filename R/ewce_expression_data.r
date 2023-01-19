#' Bootstrap cell type enrichment test for transcriptome data
#'
#' \code{ewce_expression_data} takes a differential gene expression (DGE)
#' results table and determines the probability of cell type enrichment
#' in the up- and down- regulated genes.
#'
#' @param tt Differential expression table.
#' Can be output of \link[limma]{topTable} function.
#' Minimum requirement is that one column stores a metric of
#' increased/decreased expression (i.e. log fold change, t-statistic for
#' differential expression etc) and another contains gene symbols.
#' @param sortBy Column name of metric in \code{tt}
#'  which should be used to sort up- from down- regulated genes (Default: "t").
#' @param thresh The number of up- and down- regulated genes to be included in
#' each analysis (Default: 250).
#' @param ttSpecies The species the differential expression table
#' was generated from.
#' @inheritParams bootstrap_enrichment_test
#' @inheritParams orthogene::convert_orthologs
#'
#' @return A list containing five data frames:
#' \itemize{
#'   \item \code{results}: dataframe in which each row gives the statistics
#'   (p-value, fold change and number of standard deviations from the mean)
#'   associated with the enrichment of the stated cell type in the gene list.
#'   An additional column *Direction* stores whether it the result is from the
#'   up or downregulated set.
#'   \item \code{hit.cells.up}: vector containing the summed proportion of
#'   expression in each cell type for the target list.
#'   \item \code{hit.cells.down}: vector containing the summed proportion of
#'   expression in each cell type for the target list.
#'   \item \code{bootstrap_data.up}: matrix in which each row represents the
#'   summed proportion of expression in each cell type for one of the random
#'   lists.
#'   \item \code{bootstrap_data.down}: matrix in which each row represents the
#'   summed proportion of expression in each cell type for one of the random
#'   lists.
#' }
#' @examples
#' # Load the single cell data
#' ctd <- ewceData::ctd()
#'
#' # Set the parameters for the analysis
#' # Use 3 bootstrap lists for speed, for publishable analysis use >10000
#' reps <- 3
#' # Use 5 up/down regulated genes (thresh) for speed, default is 250
#' thresh <- 5
#' annotLevel <- 1 # <- Use cell level annotations (i.e. Interneurons)
#'
#' # Load the top table
#' tt_alzh <- ewceData::tt_alzh()
#'
#' tt_results <- EWCE::ewce_expression_data(
#'     sct_data = ctd,
#'     tt = tt_alzh,
#'     annotLevel = 1,
#'     thresh = thresh,
#'     reps = reps,
#'     ttSpecies = "human",
#'     sctSpecies = "mouse"
#' )
#' @export
ewce_expression_data <- function(sct_data,
                                 annotLevel = 1,
                                 tt,
                                 sortBy = "t",
                                 thresh = 250,
                                 reps = 100,
                                 ttSpecies = NULL,
                                 sctSpecies = NULL,
                                 output_species = NULL,
                                 bg = NULL,
                                 method = "homologene",
                                 verbose = TRUE,
                                 localHub = FALSE) {
    #### Check args ####
    check_ewce_expression_data_args(sortBy = sortBy,
                                    tt = tt, 
                                    thresh = thresh) 
    #### Check species1 ###
    species <- check_species(
        genelistSpecies = output_species,
        sctSpecies = sctSpecies,
        verbose = verbose
    )
    output_species <- species$genelistSpecies
    sctSpecies <- species$sctSpecies 
    #### Check species2 ###
    species <- check_species(
        genelistSpecies = output_species,
        sctSpecies = ttSpecies,
        verbose = verbose
    )
    output_species <- species$genelistSpecies
    ttSpecies <- species$sctSpecies 
    
    #### Generate background ####  
    bg_out <- create_background_multilist(
        gene_list1 = as.character(unname(rownames(sct_data[[1]]$specificity))),
        ## Assumes 1st col contains gene names
        gene_list2 = as.character(tt[,1]),
        gene_list1_species = sctSpecies,
        gene_list2_species = ttSpecies,
        output_species = output_species,
        bg = bg,
        use_intersect = TRUE,
        method = method,
        verbose = verbose
    )
    bg <- bg_out$bg 
    #### Standardise CTD ####
    messager("Standardising sct_data.", v = verbose)
    sct_data <- standardise_ctd(
        ctd = sct_data,
        input_species = sctSpecies,
        output_species = output_species,
        force_standardise = sctSpecies!=output_species,
        dataset = "sct_data",
        method = method,
        verbose = FALSE
    )
    sctSpecies <- output_species 
    #### Convert tt orthologs ####
    tt_list <- prepare_tt(tt = tt, 
                          ttSpecies = ttSpecies, 
                          output_species = output_species, 
                          verbose = verbose)
    tt2 <- tt_list$tt; 
    tt_genecol <- tt_list$tt_genecol; 
    ttSpecies <- tt_list$ttSpecies;  
    
    #### Sort from down-->up regulated ###3
    # Sort by t-statistic
    tt3 <- tt2[order(tt2[, sortBy]), ] 
    #### Select the up/down-regulated gene sets ####
    upreg.hits <- unique(tt3[
        dim(tt3)[1]:(dim(tt3)[1] - thresh),
        tt_genecol
    ])
    downreg.hits <- unique(tt3[seq_len(thresh), tt_genecol]) 
    
    #### Run bootstrap_enrichment_test ####
    full_res_up <- bootstrap_enrichment_test(
        sct_data = sct_data,
        hits = upreg.hits,
        bg = bg,
        reps = reps,
        annotLevel = annotLevel,
        geneSizeControl = FALSE,
        genelistSpecies = ttSpecies,
        sctSpecies = sctSpecies,
        localHub = localHub
    )
    full_res_down <-
        bootstrap_enrichment_test(
            sct_data = sct_data,
            hits = downreg.hits,
            bg = bg, 
            reps = reps,
            annotLevel = annotLevel,
            geneSizeControl = FALSE,
            genelistSpecies = ttSpecies,
            sctSpecies = sctSpecies,
            localHub = localHub
        )

    joint_results <- rbind(
        cbind(full_res_up$results, Direction = "Up"),
        cbind(full_res_down$results, Direction = "Down")
    )
    hit.cells.up <- full_res_up$hit.cells
    hit.cells.down <- full_res_down$hit.cells
    bootstrap_data.up <- full_res_up$bootstrap_data
    bootstrap_data.down <- full_res_down$bootstrap_data

    return(list(
        joint_results = joint_results, 
        hit.cells.up = hit.cells.up,
        hit.cells.down = hit.cells.down,
        bootstrap_data.up = bootstrap_data.up,
        bootstrap_data.down = bootstrap_data.down
    ))
}
