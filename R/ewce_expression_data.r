#' Bootstrap cell type enrichment test for transcriptome data
#'
#' \code{ewce_expression_data} takes a differential gene expression (DGE)
#' results table and determines the probability of cell type enrichment
#' in the up- and down- regulated genes.
#'
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
#'
#' @return A list containing five data frames:
#' \itemize{
#'   \item \code{results}: dataframe in which each row gives the statistics
#'   (p-value, fold change and number of standard deviations from the mean)
#'   associated with the enrichment of the stated cell type in the gene list.
#'   An additional column *Direction* stores whether it the result is from the
#'   up or downregulated set.
#'   \item \code{hit.cells.up}: vector containing the summed proportion of
#'   expression in each cell type for the target list
#'   \item \code{hit.cells.down}: vector containing the summed proportion of
#'   expression in each cell type for the target list#'
#'   \item \code{bootstrap_data.up}: matrix in which each row represents the
#'   summed proportion of expression in each cell type for one of the random
#'   lists
#'   \item \code{bootstrap_data.down}: matrix in which each row represents the
#'   summed proportion of expression in each cell type for one of the random
#'   lists
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
    ttSpecies = "mouse",
    sctSpecies = "mouse") {
    err_msg <- paste0(
        "ERROR: tt does not contain a column with value",
        " passed in sortBy argument"
    )
    # Check the arguments
    if (!sortBy %in% colnames(tt)) {
        stop(err_msg)
    }
    err_msg2 <- paste0(
        "ERROR: length of table is less than",
        " twice the size of threshold"
    )
    if (dim(tt)[1] < (thresh * 2)) {
        stop(err_msg2)
    }

    err_msg3 <- paste0(
        "ERROR: if ttSpecies==human then there must be an ",
        "HGNC.symbol column"
    )
    err_msg4 <- paste0(
        "ERROR: if ttSpecies==human then there must be an ",
        "MGI.symbol column"
    )
    # Check that the top table has correct columns
    if (ttSpecies == "human" & !"HGNC.symbol" %in% colnames(tt)) {
        stop(err_msg3)
    }
    if (ttSpecies == "mouse" & !"MGI.symbol" %in% colnames(tt)) {
        stop(err_msg4)
    }

    if (ttSpecies == "human") {
        tt$MGI.symbol <- tt$HGNC.symbol
    }
    tt$MGI.symbol <- as.character(tt$MGI.symbol)
    tt2 <- tt

    # Sort from down-->up regulated
    tt3 <- tt2[order(tt2[, sortBy]), ] # Sort by t-statistic

    # Select the up/down-regulated gene sets
    mouse.upreg.hits <- unique(tt3[
        dim(tt3)[1]:(dim(tt3)[1] - thresh),
        "MGI.symbol"
    ])
    mouse.downreg.hits <- unique(tt3[seq_len(thresh), "MGI.symbol"])

    # Select the background gene set
    mouse.bg <- unique(tt3$MGI.symbol)

    # Do EWCE analysis
    full_res_up <- bootstrap_enrichment_test(
        sct_data = sct_data,
        hits = mouse.upreg.hits,
        bg = mouse.bg, reps = reps,
        annotLevel = annotLevel,
        geneSizeControl = FALSE,
        genelistSpecies = ttSpecies,
        sctSpecies = sctSpecies
    )
    full_res_down <-
        bootstrap_enrichment_test(
            sct_data = sct_data,
            hits = mouse.downreg.hits,
            bg = mouse.bg, reps = reps,
            annotLevel = annotLevel,
            geneSizeControl = FALSE,
            genelistSpecies = ttSpecies,
            sctSpecies = sctSpecies
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
        joint_results = joint_results, hit.cells.up = hit.cells.up,
        hit.cells.down = hit.cells.down,
        bootstrap_data.up = bootstrap_data.up,
        bootstrap_data.down = bootstrap_data.down
    ))
}
