#' Bootstrap celltype enrichment test
#'
#' \code{bootstrap.enrichment.test} takes a genelist and a single cell type 
#' transcriptome dataset and determines the probability of enrichment and fold 
#' changes for each cell type.
#'
#' @param sct_data List generated using \code{\link{generate.celltype.data}}
#' @param hits Array of MGI gene symbols containing the target gene list. 
#' Must be HGNC symbols if geneSizeControl=TRUE
#' @param bg Array of MGI gene symbols containing the background gene list. 
#' Must be HGNC symbols if geneSizeControl=TRUE
#' @param genelistSpecies Either 'mouse' or 'human' depending on whether MGI 
#' or HGNC symbols are used for gene lists
#' @param sctSpecies  Either 'mouse' or 'human' depending on whether MGI or 
#' HGNC symbols are used for the single cell dataset
#' @param reps Number of random gene lists to generate (default=100 but should 
#' be over 10000 for publication quality results)
#' @param annotLevel an integer indicating which level of the annotation to 
#' analyse. Default = 1.
#' @param geneSizeControl a logical indicating whether you want to control for 
#' GC content and transcript length. Recommended if the gene list originates 
#' from genetic studies. Default is FALSE. If set to TRUE then human gene lists 
#' should be used rather than mouse.
#' @param controlledCT (optional) If not NULL, and instead is the name of a 
#' cell type, then the bootstrapping controls for expression within that 
#' cell type
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
#' @examples
#' library(ewceData)
#' # Load the single cell data
#' ctd <- ctd()
#'
#' # Set the parameters for the analysis
#' # Use 100 bootstrap lists for speed, for publishable analysis use >10000
#' reps <- 100 
#'
#' # Load the gene list and get human orthologs
#' example_genelist <- example_genelist()
#' mouse_to_human_homologs <- mouse_to_human_homologs()
#' m2h <- unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])
#' mouse.hits <- 
#'     unique(m2h[m2h$HGNC.symbol %in% example_genelist, "MGI.symbol"])
#' human.hits <- 
#'     unique(m2h[m2h$HGNC.symbol %in% example_genelist, "HGNC.symbol"])
#' human.bg <- unique(m2h$HGNC.symbol)
#' mouse.bg <- unique(m2h$MGI.symbol)
#'
#' # Bootstrap significance test, no control for transcript length or GC content
#' full_results <- bootstrap.enrichment.test(
#'     sct_data = ctd, hits = mouse.hits,
#'     bg = mouse.bg, reps = reps, annotLevel = 2, sctSpecies = "mouse", 
#'     genelistSpecies = "mouse"
#' )
#'
#' # Bootstrap significance test control for transcript length and GC content
#' full_results <- bootstrap.enrichment.test(
#'     sct_data = ctd, hits = human.hits,
#'     bg = human.bg, reps = reps, annotLevel = 2, geneSizeControl = TRUE,
#'     sctSpecies = "mouse", genelistSpecies = "human"
#' )
#' @export
#' @import stats
# @importFrom reshape2 melt
# @import plyr
bootstrap.enrichment.test <- function(sct_data = NA, hits = NA, bg = NA,
                                        genelistSpecies = "mouse",
                                        sctSpecies = "mouse", reps = 100,
                                        annotLevel = 1, geneSizeControl = FALSE,
                                        controlledCT = NULL) {
    checkedLists <- check.ewce.genelist.inputs(
        sct_data, hits, bg, genelistSpecies,
        sctSpecies, geneSizeControl
    )
    hits <- checkedLists$hits
    bg <- checkedLists$bg

    # Check an SCT dataset was provided
    if (unique(is.na(sct_data))) {
        stop("ERROR: must provide valid single cell dataset")
    }
    err_msg <- paste0("ERROR: invalid celltype name passed in controlledCT.", 
                        " This argument is optional. Leave empty if you do not",
                        " wish to control for a celltypes expression.")
    # Check if controlling for another celltype
    if (!is.null(controlledCT)) {
        if (!controlledCT %in% colnames(sct_data[[1]]$specificity)) {
            stop(err_msg)
        }
    }
    err_msg2 <- paste0("ERROR! AFTER CALCULATING BOOTSTRAPPING NETWORK",
                        " WITH LENGTH + GC CONTROLS, size of list_network",
                        " is not same length as hitGenes")
    # IF USING GENESIZE AND GC-CONTENT MATCHING, THEN GENERATE THE SAMPLE LISTS
    if (geneSizeControl == TRUE) {
        control_related <- prepare.genesize.control.network(
            hits = hits,
            bg = bg,
            numBOOT = reps,
            sctSpecies = sctSpecies
        )
        control_network <- control_related[["list_network"]]
        hitGenes <- control_related[["hitGenes"]]
        nonHits <- unique(control_related[["list_network"]]) # mouse.bg
        combinedGenes <- c(hitGenes, nonHits) # c(mouse.hits,mouse.bg)
        if (length(hitGenes) != dim(control_network)[2]) {
            stop(err_msg2)
        }
    } else {
        hitGenes <- hits # mouse.hits
        nonHits <- bg # mouse.bg
        combinedGenes <- c(hits, bg) # c(mouse.hits,mouse.bg)
    }


    if (!is.null(names(sct_data))) {
        # This is necessary in case further meta-data such as $name is used
        numLevels <- sum(names(sct_data) == "") 
    } else {
        numLevels <- length(sct_data)
    }
    for (lv in seq_len(numLevels)) {
        gN <- rownames(sct_data[[lv]]$mean_exp)
        sct_data[[lv]]$mean_exp <- 
            sct_data[[lv]]$mean_exp[gN %in% combinedGenes, ]
        sct_data[[lv]]$specificity <- 
            sct_data[[lv]]$specificity[gN %in% combinedGenes, ]
    }

    # GET WEIGHTING FOR EACH CELL TYPE FOR HIT LIST (hit.cells) AND A MATRIX TO 
    #STORE THE BOOTSTRAP WEIGHTINGS (bootstrap_data)
    cells <- unique(colnames(sct_data[[annotLevel]]$specificity))

    # GENERATE hit.cells and bootstrap_data IN ONE GO
    print(hitGenes)
    if (!geneSizeControl) {
        control_network <- NULL
    }
    sumProp <- get_summed_proportions(
        hitGenes,
        sct_data,
        annotLevel,
        reps,
        geneSizeControl,
        controlledCT,
        control_network = control_network
    )

    hit.cells <- sumProp$hit.cells
    bootstrap_data <- sumProp$bootstrap_data

    # EXTRACTING THE DETAILS:
    # - CALCULATING P-VALUE, FOLD CHANGE AND MARKERS ETC
    count <- 0
    for (ct in cells) {
        print(ct)
        count <- count + 1
        # For cell type 'ct' get the bootstrap and target list values
        ct_boot_dist <- 
            bootstrap_data[, colnames(sct_data[[annotLevel]]$specificity) == ct]
        hit_sum <- hit.cells[colnames(sct_data[[annotLevel]]$specificity) == ct]
        # Get propability and fold change of enrichment
        p <- sum(ct_boot_dist >= hit_sum) / reps
        print(p)
        fold_change <- hit_sum / mean(ct_boot_dist)
        sd_from_mean <- (hit_sum - mean(ct_boot_dist)) / sd(ct_boot_dist)
        ct_root <- ct
        if (p < 0.05) {
            # If cell type is significant, print the contributing genes:
            print(sprintf("Fold enrichment: %s", fold_change))
            print(sprintf("Standard deviations from mean: %s", sd_from_mean))
        }
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
        print("")
    }

    full_results <- list(
        results = results,
        hit.cells = hit.cells,
        bootstrap_data = bootstrap_data
    )
    return(full_results)
}
