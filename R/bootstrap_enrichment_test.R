#' Bootstrap cell type enrichment test
#'
#' \code{bootstrap_enrichment_test} takes a genelist and a single cell type
#' transcriptome dataset and determines the probability of enrichment and fold
#' changes for each cell type.
#'
#' @param sct_data List generated using \link[EWCE]{generate_celltype_data}.
#' @param hits List of gene symbols containing the target gene list.
#' Will automatically be converted to human gene symbols
#' if \code{geneSizeControl=TRUE}.
#' @param bg List of gene symbols containing the background gene list
#' (including hit genes). If \code{bg=NULL},
#'  an appropriate gene background will be created automatically.
#' @param genelistSpecies Species that \code{hits} genes came from
#' (no longer limited to just "mouse" and "human").
#' See \link[EWCE]{list_species} for all available species.
#' @param sctSpecies  Species that \code{sct_data} came from
#' (no longer limited to just "mouse" and "human").
#' See \link[EWCE]{list_species} for all available species.
#' @param output_species Species to convert \code{sct_data} and \code{hits} to
#' (Default: "human").
#' See \link[EWCE]{list_species} for all available species.
#' @param reps Number of random gene lists to generate (\emph{Default: 100},
#'  but should be >=10,000 for publication-quality results).
#' @param no_cores Number of cores to parallelise
#' bootstrapping \code{reps} over.
#' @param annotLevel An integer indicating which level of \code{sct_data} to
#' analyse (\emph{Default: 1}).
#' @param geneSizeControl Whether you want to control for
#' GC content and transcript length. Recommended if the gene list originates
#' from genetic studies (\emph{Default: FALSE}).
#' If set to \code{TRUE}, then \code{hits} must be from humans.
#' @param controlledCT [Optional] If not NULL, and instead is the name of a
#' cell type, then the bootstrapping controls for expression within that
#' cell type.
#' @param sort_results Sort enrichment results from
#'  smallest to largest p-values.
#' @param mtc_method Multiple-testing correction method
#' (passed to \link[stats]{p.adjust}).
#' @param verbose Print messages.
#' @inheritParams orthogene::convert_orthologs
#'
#' @return A list containing three data frames:
#' \itemize{
#'   \item \code{results}: dataframe in which each row gives the statistics
#'   (p-value, fold change and number of standard deviations from the mean)
#'   associated with the enrichment of the stated cell type in the gene list
#'   \item \code{hit.cells}: vector containing the summed proportion of
#'   expression in each cell type for the target list
#'   \item \code{bootstrap_data}: matrix in which each row represents the
#'   summed proportion of expression in each cell type for one of the
#'   random lists
#' }
#'
#' @examples
#' # Load the single cell data
#' ctd <- ewceData::ctd()
#' # Set the parameters for the analysis
#' # Use 3 bootstrap lists for speed, for publishable analysis use >=10,000
#' reps <- 3
#' # Load gene list from Alzheimer's disease GWAS
#' example_genelist <- ewceData::example_genelist()
#'
#' # Bootstrap significance test, no control for transcript length or GC content
#' full_results <- EWCE::bootstrap_enrichment_test(
#'     sct_data = ctd,
#'     hits = example_genelist,
#'     reps = reps,
#'     annotLevel = 1,
#'     sctSpecies = "mouse",
#'     genelistSpecies = "human"
#' )
#' @export
#' @importFrom dplyr %>% arrange desc
#' @importFrom stats p.adjust sd
#' @importFrom orthogene create_background
bootstrap_enrichment_test <- function(sct_data = NULL,
                                      hits = NULL,
                                      bg = NULL,
                                      genelistSpecies = NULL,
                                      sctSpecies = NULL,
                                      output_species = "human",
                                      method = "homologene",
                                      reps = 100,
                                      no_cores = 1,
                                      annotLevel = 1,
                                      geneSizeControl = FALSE,
                                      controlledCT = NULL,
                                      mtc_method = "BH",
                                      sort_results = TRUE,
                                      verbose = TRUE) {
    core_allocation <- assign_cores(
        worker_cores = no_cores,
        verbose = verbose
    )
    #### Check controlledCT ####
    controlledCT <- fix_celltype_names(celltypes = controlledCT)
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
        controlledCT = controlledCT
    )
    #### Create background if none provided ####
    bg <- orthogene::create_background(
        species1 = sctSpecies,
        species2 = genelistSpecies,
        output_species = output_species,
        bg = bg,
        method = method,
        verbose = verbose
    )
    #### Convert CTD to standardized human genes ####
    messager("Standardising CellTypeDataset", v = verbose)
    sct_data <- standardise_ctd(
        ctd = sct_data,
        input_species = sctSpecies,
        output_species = output_species,
        dataset = "sct_data",
        method = method,
        verbose = FALSE
    )
    sctSpecies <- output_species
    #### Convert gene list inputs to standardized human genes ####
    checkedLists <- check_ewce_genelist_inputs(
        sct_data = sct_data,
        hits = hits,
        bg = bg,
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies,
        geneSizeControl = geneSizeControl,
        output_species = output_species,
        verbose = verbose
    )
    #### check_ewce_genelist_inputs converts hits to "human" by default ####
    genelistSpecies <- checkedLists$output_species
    hits <- checkedLists$hits
    bg <- checkedLists$bg
    #### Correct for genesize and gc-content matching ####
    # If using genesize and GC-content matching, generate the sample lists
    if (isTRUE(geneSizeControl)) {
        messager("Running with gene size control.", v = verbose)
        control_related <- prepare_genesize_control_network(
            hits = hits,
            bg = bg,
            reps = reps,
            sctSpecies = sctSpecies
        )
        control_network <- control_related[["list_network"]]
        hitGenes <- control_related[["hitGenes"]]
        nonHits <- unique(control_related[["list_network"]]) # mouse.bg
        combinedGenes <- c(hitGenes, nonHits) # c(mouse.hits,mouse.bg)
        if (length(hitGenes) != dim(control_network)[2]) {
            err_msg2 <- paste0(
                "ERROR! AFTER CALCULATING BOOTSTRAPPING NETWORK",
                " WITH LENGTH + GC CONTROLS, size of list_network",
                " is not same length as hitGenes"
            )
            stop(err_msg2)
        }
    } else {
        messager("Running without gene size control.", v = verbose)
        hitGenes <- hits # mouse.hits
        nonHits <- bg # mouse.bg
        combinedGenes <- c(hits, nonHits) # c(mouse.hits,mouse.bg)
    }

    lvls <- get_ctd_levels(ctd = sct_data)
    numLevels <- length(lvls)
    for (lv in seq_len(numLevels)) {
        gN <- rownames(sct_data[[lv]]$mean_exp)
        sct_data[[lv]]$mean_exp <-
            sct_data[[lv]]$mean_exp[gN %in% combinedGenes, ]
        sct_data[[lv]]$specificity <-
            sct_data[[lv]]$specificity[gN %in% combinedGenes, ]
    }

    # GET WEIGHTING FOR EACH CELL TYPE FOR HIT LIST (hit.cells) AND A MATRIX TO
    # STORE THE BOOTSTRAP WEIGHTINGS (bootstrap_data)
    cells <- unique(colnames(sct_data[[annotLevel]]$specificity))
    # GENERATE hit.cells and bootstrap_data IN ONE GO
    messager(formatC(length(unique(hitGenes)), big.mark = ","),
        "hit genes remain after filtering.",
        v = verbose
    )
    if (!geneSizeControl) {
        control_network <- NULL
    }
    sumProp <- get_summed_proportions(
        hitGenes = hitGenes,
        sct_data = sct_data,
        annotLevel = annotLevel,
        reps = reps,
        no_cores = no_cores,
        geneSizeControl = geneSizeControl,
        controlledCT = controlledCT,
        control_network = control_network,
        verbose = verbose
    )
    hit.cells <- sumProp$hit.cells
    bootstrap_data <- sumProp$bootstrap_data

    #### EXTRACTING THE DETAILS ####
    # - CALCULATING P-VALUE, FOLD CHANGE AND MARKERS ETC
    messager("Testing for enrichment in", length(cells), "cell types...",
        v = verbose
    )
    count <- 0
    for (ct in cells) {
        count <- count + 1
        # For cell type 'ct' get the bootstrap and target list values
        ct_boot_dist <-
            bootstrap_data[, colnames(sct_data[[annotLevel]]$specificity) == ct]
        hit_sum <- hit.cells[colnames(sct_data[[annotLevel]]$specificity) == ct]
        # Get probability and fold change of enrichment
        p <- sum(ct_boot_dist >= hit_sum) / reps
        # messager(p,v=verbose)
        fold_change <- hit_sum / mean(ct_boot_dist)
        sd_from_mean <- (hit_sum - mean(ct_boot_dist)) /
            stats::sd(ct_boot_dist)
        ct_root <- ct
        if (count == 1) {
            results <- data.frame(
                CellType = ct,
                annotLevel = annotLevel,
                p = p,
                fold_change = fold_change,
                sd_from_mean = sd_from_mean
            )
        } else {
            results <- rbind(
                results,
                data.frame(
                    CellType = ct,
                    annotLevel = annotLevel,
                    p = p,
                    fold_change = fold_change,
                    sd_from_mean = sd_from_mean
                )
            )
        }
    } ## End for loop

    #### Sort results ####
    if (sort_results) {
        messager("Sorting results by p-value.", v = verbose)
        results <- dplyr::arrange(results, p, dplyr::desc(sd_from_mean))
    }
    #### Apply multiple testing correction ####
    if (!is.null(mtc_method)) {
        messager("Computing", paste0(mtc_method, "-corrected"), "q-values.",
            v = verbose
        )
        results$q <- stats::p.adjust(
            p = results$p,
            method = mtc_method
        )
    }
    #### Report results ####
    report_results(
        results = results,
        verbose = verbose
    )
    #### Return results list ####
    full_results <- list(
        results = results,
        hit.cells = hit.cells,
        bootstrap_data = bootstrap_data
    )
    return(full_results)
}
