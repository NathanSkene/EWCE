#' Generate bootstrap plots
#'
#' \code{generate_bootstrap_plots} takes a gene list and a single cell type
#' transcriptome dataset and generates plots which show how the expression of
#' the genes in the list compares to those in randomly generated gene lists.
#' @param full_results The full output of
#' \link[EWCE]{bootstrap_enrichment_test} for the same gene list.
#' @param listFileName String used as the root for files saved using this
#' function.
#' @param save_dir Directory where the BootstrapPlots folder should be saved,
#' default is a temp directory.
#' @param adj_pval_thresh Adjusted p-value threshold of celltypes to include
#' in plots.
#' @inheritParams bootstrap_enrichment_test
#' @inheritParams bootstrap_plot
#' @inheritParams orthogene::create_background
#' @inheritParams ggplot2::facet_grid
#'
#' @returns Saves a set of pdf files containing graphs and returns the file where
#' they are saved. These will be saved with the file name adjusted using the
#' value of \code{listFileName}. The files are saved into the 
#' 'BootstrapPlot' folder.
#' Files start with one of the following:
#' \itemize{
#'   \item \code{qqplot_noText}: sorts the gene list according to how enriched
#'   it is in the relevant cell type. Plots the value in the target list against
#'   the mean value in the bootstrapped lists.
#'   \item \code{qqplot_wtGSym}: as above but labels the gene symbols for the
#'   highest expressed genes.
#'   \item \code{bootDists}: rather than just showing the mean of the
#'   bootstrapped lists, a boxplot shows the distribution of values
#'   \item \code{bootDists_LOG}: shows the bootstrapped distributions with the
#'   y-axis shown on a log scale
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom orthogene create_background
#' @examples
#' ## Load the single cell data
#' sct_data <- ewceData::ctd()
#'
#' ## Set the parameters for the analysis
#' ## Use 5 bootstrap lists for speed, for publishable analysis use >10000
#' reps <- 5
#'
#' ## Load the gene list and get human orthologs
#' hits <- ewceData::example_genelist()
#'
#' ## Bootstrap significance test,
#' ##  no control for transcript length or GC content
#' ## Use pre-computed results to speed up example
#' full_results <- EWCE::example_bootstrap_results()
#'
#' ### Skip this for example purposes
#' # full_results <- EWCE::bootstrap_enrichment_test(
#' #    sct_data = sct_data,
#' #    hits = hits,
#' #    reps = reps,
#' #    annotLevel = 1,
#' #    sctSpecies = "mouse",
#' #    genelistSpecies = "human"
#' # )
#'
#' output <- EWCE::generate_bootstrap_plots(
#'     sct_data = sct_data,
#'     hits = hits,
#'     reps = reps,
#'     full_results = full_results,
#'     sctSpecies = "mouse",
#'     genelistSpecies = "human",
#'     annotLevel = 1
#' )
generate_bootstrap_plots <- function(sct_data = NULL,
                                     hits = NULL,
                                     bg = NULL,
                                     genelistSpecies = NULL,
                                     sctSpecies = NULL,
                                     output_species = "human",
                                     method = "homologene",
                                     reps = 100,
                                     annotLevel = 1,
                                     geneSizeControl = FALSE,
                                     full_results = NULL,
                                     listFileName = paste0("_level",
                                                           annotLevel),
                                     adj_pval_thresh = 0.05,
                                     facets = "CellType",
                                     scales = "free_x",
                                     save_dir = file.path(tempdir(),
                                                          "BootstrapPlots"),  
                                     show_plot = TRUE,
                                     # min_genes = 4,
                                     verbose = TRUE) {
    # devoptera::args2vars(generate_bootstrap_plots)
    # #' @param min_genes The minimum number of genes in \code{hits} that are 
    # #' also in the single cell dataset & background gene set.
    #### Set min_genes ####
    min_genes <- Sys.getenv("min_genes")
    min_genes <- if(min_genes=="") 4 else as.numeric(min_genes) 
    
    #### Check species ####
    species <- check_species(
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies
    )
    genelistSpecies <- species$genelistSpecies
    sctSpecies <- species$sctSpecies
    sctSpecies_origin <- species$sctSpecies_origin
    #### Check bootstrap args ####
    check_bootstrap_args(
        sct_data = sct_data,
        hits = hits,
        annotLevel = annotLevel,
        reps = reps,
        fix_celltypes = TRUE
    )
    
    #### Create background if none provided ####
    #if statement added down to issue with orthogene:
    #https://github.com/neurogenomics/orthogene/issues/22
    if(is.null(bg) | 
       !all(list(sctSpecies,
                 genelistSpecies,
                 sctSpecies_origin)==output_species)){
      bg <- orthogene::create_background(
        species1 = sctSpecies_origin,
        species2 = genelistSpecies,
        output_species = output_species,
        method = method,
        bg = bg,
        verbose = verbose
      )
    }else{
      bg <- unique(bg)
    }
    #### Standardise CTD ####
    messager("Standardising sct_data.", v = verbose)
    sct_data <- standardise_ctd(
        ctd = sct_data,
        input_species = sctSpecies,
        output_species = output_species,
        dataset = "sct_data",
        method = method,
        verbose = FALSE
    )
    sctSpecies <- output_species
    #### Check full results AFTER ctd has been standardised ####
    full_results <- fix_celltype_names_full_results(
        full_results = full_results,
        verbose = verbose
    )
    check_full_results(
        full_results = full_results,
        sct_data = sct_data
    ) 
    results <- full_results$results
    #### Check gene lists ####
    checkedLists <- check_ewce_genelist_inputs(
        sct_data = sct_data,
        hits = hits,
        bg = bg,
        sctSpecies = sctSpecies,
        genelistSpecies = genelistSpecies,
        sctSpecies_origin = sctSpecies_origin,
        geneSizeControl = geneSizeControl,
        output_species = output_species,
        min_genes = min_genes,
        verbose = verbose
    )
    hits <- checkedLists$hits
    bg <- checkedLists$bg
    combinedGenes <- unique(c(hits, bg))
    #### Check significant results #### 
    signif_ct <- rownames(results)[results$q < adj_pval_thresh]
    if(length(signif_ct)==0){
        stp <- "Must have >0 significant celltypes."
        stop(stp)
    } else {
        messager(length(signif_ct),"celltype(s) remain @ <=",adj_pval_thresh,
                 v=verbose)  
    } 
    #### Get specificity data of bootstrapped genes ####
    ## Still need to regenerate these exp_mats for fig4, 
    ## since this info is not stored in the full_results.
    exp_mats <- generate_bootstrap_plots_exp_mats(sct_data=sct_data,
                                                  annotLevel=annotLevel, 
                                                  reps=reps,
                                                  combinedGenes=combinedGenes,
                                                  hits=hits,
                                                  verbose=verbose)
    #### Use precomputed gene_data if available ####
    if(!is.null(full_results) &&
       all(!is.na(full_results)) &&
       !is.null(full_results$gene_data)){
        gene_data <- full_results$gene_data
    } else { 
        cgs <- compute_gene_scores(sct_data = sct_data, 
                                   annotLevel = annotLevel,  
                                   hits = hits, 
                                   reps = reps,
                                   combinedGenes = combinedGenes,
                                   exp_mats = exp_mats,
                                   return_hit_exp = TRUE,
                                   verbose = verbose)  
        gene_data <- cgs$gene_data
    } 
    #### Iteratively create QQ plots #### 
    messager("Generating bootstrap plot for",
             length(signif_ct),"celltype(s).", v = verbose)
    out <- bootstrap_plot(
        gene_data = gene_data, 
        exp_mats = exp_mats,
        save_dir = save_dir,
        listFileName = listFileName, 
        facets = facets,
        scales = scales,
        show_plot = show_plot,
        verbose = verbose
    ) 
    return(out)
}
