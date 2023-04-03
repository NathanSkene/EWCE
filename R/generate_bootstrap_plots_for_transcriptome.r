#' Generate bootstrap plots
#'
#' Takes a gene list and a single cell type transcriptome dataset 
#' and generates plots which show how the expression of the genes in
#'  the list compares to those in randomly generated gene lists.
#'
#' @param full_results The full output of
#' \link[EWCE]{ewce_expression_data} for the same gene list.
#' @param listFileName String used as the root for files saved using
#' this function. 
#' @param showGNameThresh Integer. If a gene has over X percent of it's
#' expression proportion in a cell type, then list the gene name.
#' @param sig_only Should plots only be generated for cells which have
#' significant changes?
#' @param sig_col Column name in \code{tt} that contains the
#'  significance values.
#' @param sig_thresh Threshold by which to filter \code{tt} by \code{sig_col}.
#' @param celltype_col Column within \code{tt} that contains celltype names.
#' @param plot_types Plot types to generate. 
#' @param save_dir Directory where the BootstrapPlots folder should be saved,
#' default is a temp directory.
#' @inheritParams bootstrap_enrichment_test
#' @inheritParams ewce_expression_data
#' @inheritParams orthogene::convert_orthologs
#'
#' @returns Saves a set of PDF files containing graphs. 
#' Then returns a nested list with each \code{plot} and
#'  the \code{path} where it was saved to.
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
#' @examples
#' ## Load the single cell data
#' ctd <- ewceData::ctd()
#'
#' ## Set the parameters for the analysis
#' ## Use 3 bootstrap lists for speed, for publishable analysis use >10,000
#' reps <- 3
#' annotLevel <- 1 # <- Use cell level annotations (i.e. Interneurons)
#' ## Use 5 up/down regulated genes (thresh) for speed, default is 250
#' thresh <- 5
#'
#' ## Load the top table
#' tt_alzh <- ewceData::tt_alzh()
#'
#' ## See ?example_transcriptome_results for full code to produce tt_results
#' tt_results <- EWCE::example_transcriptome_results()
#'
#' ## Bootstrap significance test,
#' ## no control for transcript length or GC content
#' savePath <- EWCE::generate_bootstrap_plots_for_transcriptome(
#'     sct_data = ctd,
#'     tt = tt_alzh,
#'     thresh = thresh,
#'     annotLevel = 1,
#'     full_results = tt_results,
#'     listFileName = "examples",
#'     reps = reps,
#'     ttSpecies = "human",
#'     sctSpecies = "mouse", 
#'     # Only do one plot type for demo purposes
#'     plot_types = "bootstrap" 
#' ) 
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
generate_bootstrap_plots_for_transcriptome <- function(
    sct_data,
    tt, 
    bg = NULL,
    thresh = 250,
    annotLevel = 1,
    reps = 100,
    full_results = NA,
    listFileName = "",
    showGNameThresh = 25,
    ttSpecies = NULL,
    sctSpecies = NULL,  
    output_species = NULL,
    sortBy = "t",
    sig_only = TRUE,
    sig_col = "q",
    sig_thresh = 0.05,
    celltype_col = "CellType",
    plot_types = c("bootstrap",
                   "bootstrap_distributions",
                   "log_bootstrap_distributions"),
    save_dir = file.path(tempdir(),"BootstrapPlots"),
    method = "homologene",
    verbose = TRUE) { 
    #### Check inputs ####
    plot_types <- tolower(plot_types)
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
    #### Fix celltype names ####
    full_results <- fix_celltype_names_full_results(full_results = full_results)
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
    # sct_genes <- unname(bg_out$gene_list1)
    # tt_genes <- unname(bg_out$gene_list2)
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
    #### Check args ####
    check_args_for_bootstrap_plot_generation(
        sct_data = sct_data,
        tt = tt,
        thresh = thresh,
        annotLevel = annotLevel,
        reps = reps,
        full_results = full_results,
        listFileName = listFileName,
        showGNameThresh = showGNameThresh, 
        sortBy = sortBy
    )
    #### Convert tt orthologs ####
    tt_list <- prepare_tt(tt = tt, 
                          ttSpecies = ttSpecies, 
                          output_species = output_species, 
                          verbose = verbose)
    tt <- tt_list$tt; 
    tt_genecol <- tt_list$tt_genecol; 
    ttSpecies <- tt_list$ttSpecies;  
    
    ### Create plots of up/down regulated genes in each celltype ####
    for (dirS in c("Up", "Down")) {
        #### Sort tt by up/down regulated ####
        a <- full_results$joint_results
        results <- a[as.character(a$Direction) == dirS, ]  
        if (dirS == "Up") {
            tt <- tt[order(tt[, sortBy], decreasing = TRUE), ]
        }
        if (dirS == "Down") {
            tt <- tt[order(tt[, sortBy], decreasing = FALSE), ]
        }
        #### Drop hits genes not in expression data #### 
        ### IMPORTANT!: Keep in the this specific order 
        {
            #### sct genes ####
            spec <- sct_data[[annotLevel]]$specificity
            sct_genes <- unique(as.character(unname(rownames(spec))))
            #### hits ####
            hits <- unique( as.character(unname(tt[,tt_genecol])) )
            hits <- hits[hits %in% sct_genes]
            bg <- bg[!bg %in% hits] 
            bg <- bg[bg %in% sct_genes]
            #### Combined genes ####
            combinedGenes <- unique(c(hits, bg))
            combinedGenes <- unique(as.character(unname(combinedGenes)))
        }
        #### Get expression data of bootstrapped genes ####
        if (isTRUE(sig_only)) {
            signif_res <- as.character(results[[celltype_col]])[
                results[[sig_col]] < sig_thresh]
        } else {
            signif_res <- as.character(results[[celltype_col]])
        }
        signif_res <- fix_celltype_names(celltypes = signif_res)
        #### Create matrices of bootstrapped genes ####
        exp_mats <- get_exp_data_for_bootstrapped_genes(
            results = results,
            signif_res = signif_res,
            sct_data = sct_data,
            hits  = hits,
            combinedGenes = combinedGenes,
            annotLevel = annotLevel,
            nReps = reps
        ) 
        #### Get expression levels of the hit genes ####
        hit.exp <- sct_data[[annotLevel]]$specificity[hits, ]
        graph_theme <- theme_graph() 
        tag <- sprintf("thresh%s__dir%s", thresh, dirS)
        
        plots <- list()
        #### Create QQ plots ####
        for (cc in signif_res) {
            mean_boot_exp <- apply(exp_mats[[cc]], 2, mean)
            hit_exp <- sort(hit.exp[, cc])
            hit_exp_names <- rownames(hit.exp)[order(hit.exp[, cc])]
            dat <- data.frame(
                boot = mean_boot_exp,
                hit = hit_exp,
                Gnames = hit_exp_names
            )
            dat$hit <- dat$hit * 100
            dat$boot <- dat$boot * 100 
            maxX <- max(dat$boot, na.rm = TRUE) +
                0.1 * max(dat$boot, na.rm = TRUE)
            #### Plot several variants of the graph #### 
            if("bootstrap" %in% plot_types){
                plots[[cc]][["bootstrap"]] <- 
                    bootstrap_plots_for_transcriptome(
                    dat = dat,
                    tag = tag,
                    listFileName = listFileName,
                    cc = cc,
                    showGNameThresh = showGNameThresh,
                    graph_theme = graph_theme,
                    maxX = maxX,
                    save_dir = save_dir
                )
            } 
            #### Plot with bootstrap distribution ####
            if("bootstrap_distributions" %in% plot_types){
                plots[[cc]][["bootstrap_distributions"]] <- 
                    plot_with_bootstrap_distributions(
                    exp_mats = exp_mats,
                    cc = cc,
                    hit_exp = hit_exp,
                    tag = tag,
                    listFileName = listFileName,
                    graph_theme = graph_theme,
                    save_dir = save_dir
                )
            }
            #### Plot with LOG bootstrap distribution ####
            if("log_bootstrap_distributions" %in% plot_types){
                plots[[cc]][["log_bootstrap_distributions"]] <- 
                    plot_log_bootstrap_distributions(
                    dat = dat,
                    exp_mats = exp_mats,
                    cc = cc,
                    hit_exp = hit_exp,
                    tag = tag,
                    listFileName = listFileName,
                    graph_theme = graph_theme,
                    save_dir = save_dir
                )
            }
        }
    }
    #### Return nested list of plots and paths ####
    return(plots)
}
