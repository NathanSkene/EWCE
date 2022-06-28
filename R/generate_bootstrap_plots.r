#' Generate bootstrap plots
#'
#' \code{generate_bootstrap_plots} takes a gene list and a single cell type
#' transcriptome dataset and generates plots which show how the expression of
#' the genes in the list compares to those in randomly generated gene lists
#'
#'
#' @param full_results The full output of
#' \link[EWCE]{bootstrap_enrichment_test} for the same gene list.
#' @param listFileName String used as the root for files saved using this
#' function.
#' @param savePath Directory where the BootstrapPlots folder should be saved,
#' default is a temp directory.
#' @inheritParams bootstrap_enrichment_test
#' @inheritParams orthogene::create_background
#'
#' @return Saves a set of pdf files containing graphs and returns the file where
#' they are saved. These will be saved with the filename adjusted using the
#' value of listFileName. The files are saved into the 'BootstrapPlot' folder.
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
#'
#' @examples
#' ## Load the single cell data
#' ctd <- ewceData::ctd()
#'
#' ## Set the parameters for the analysis
#' ## Use 5 bootstrap lists for speed, for publishable analysis use >10000
#' reps <- 5
#'
#' ## Load the gene list and get human orthologs
#' hits <- ewceData::example_genelist()[1:100]
#'
#' ## Bootstrap significance test,
#' ##  no control for transcript length or GC content
#' ## Use pre-computed results to speed up example
#' full_results <- EWCE::example_bootstrap_results()
#'
#' ### Skip this for example purposes
#' # full_results <- EWCE::bootstrap_enrichment_test(
#' #    sct_data = ctd,
#' #    hits = hits,
#' #    reps = reps,
#' #    annotLevel = 1,
#' #    sctSpecies = "mouse",
#' #    genelistSpecies = "human"
#' # )
#'
#' plot_file_path <- EWCE::generate_bootstrap_plots(
#'     sct_data = ctd,
#'     hits = hits,
#'     reps = reps,
#'     full_results = full_results,
#'     listFileName = "Example",
#'     sctSpecies = "mouse",
#'     genelistSpecies = "human",
#'     annotLevel = 1,
#'     savePath = tempdir()
#' )
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom orthogene create_background
generate_bootstrap_plots <- function(sct_data = NULL,
                                     hits = NULL,
                                     bg = NULL,
                                     genelistSpecies = NULL,
                                     sctSpecies = NULL,
                                     output_species = "human",
                                     method = "homologene",
                                     reps = 100,
                                     annotLevel = 1,
                                     full_results = NA,
                                     listFileName = "",
                                     savePath = tempdir(),
                                     verbose = TRUE) {
    #### Check species ####
    species <- check_species(
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies
    )
    genelistSpecies <- species$genelistSpecies
    sctSpecies <- species$sctSpecies
    #### Check bootstrap args ####
    check_bootstrap_args(
        sct_data = sct_data,
        hits = hits,
        annotLevel = annotLevel,
        reps = reps,
        fix_celltypes = TRUE
    )
    #### Create background if none provided ####
    bg <- orthogene::create_background(
        species1 = sctSpecies,
        species2 = genelistSpecies,
        output_species = output_species,
        bg = bg,
        verbose = verbose
    )
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
    #### Add annotLevel to file name tag ###
    listFileName <- sprintf("%s_level%s", listFileName, annotLevel)
    results <- full_results$results
    #### Check gene lists ####
    checkedLists <- check_ewce_genelist_inputs(
        sct_data = sct_data,
        hits = hits,
        bg = bg,
        sctSpecies = sctSpecies,
        genelistSpecies = genelistSpecies,
        verbose = FALSE
    )
    hits <- checkedLists$hits
    bg <- checkedLists$bg
    combinedGenes <- unique(c(hits, bg))
    #### Get specificity data of bootstrapped genes ####
    signif_res <- rownames(results)[results$p < 0.05]
    nReps <- reps
    exp_mats <- list()
    for (cc in signif_res) {
        exp_mats[[cc]] <- matrix(0,
            nrow = nReps,
            ncol = length(hits)
        )
        rownames(exp_mats[[cc]]) <- sprintf("Rep%s", seq_len(nReps))
    }
    for (s in seq_len(nReps)) {
        bootstrap_set <- sample(x = combinedGenes, 
                                size = length(hits))
        ValidGenes <- rownames(sct_data[[annotLevel]]$specificity)[
            rownames(sct_data[[annotLevel]]$specificity) %in% bootstrap_set
        ]
        expD <- sct_data[[annotLevel]]$specificity[ValidGenes, ]
        for (cc in signif_res) {
            exp_mats[[cc]][s, ] <- sort(expD[, cc])
        }
    }
    #### Get specificity scores of the hit genes ####
    hit.exp <- sct_data[[annotLevel]]$specificity[hits, ]
    #### Create subdir ####
    boot_dir <- file.path(savePath,"BootstrapPlots")
    if (!file.exists(boot_dir)) {
        dir.create(file.path(savePath, "BootstrapPlots"))
    }
    #### Iteratively create QQ plots ####
    for (cc in signif_res) {
        messager("Generating bootstrap plot:", cc, v = verbose)
        bootstrap_plot(
            exp_mats = exp_mats,
            hit.exp = hit.exp,
            cc = cc,
            savePath = savePath,
            listFileName = listFileName
        )
    }
    #### return path to the saved directory in case tempdir() used ####
    return(savePath)
}
