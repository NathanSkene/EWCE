#' Multiple EWCE results from multiple studies
#'
#' \code{merged_ewce} combines enrichment results from multiple studies 
#' targetting the same scientific problem
#'
#' @param results a list of EWCE results generated using 
#' \code{\link{add.res.to.merging.list}}
#' @param reps Number of random gene lists to generate (default=100 but should 
#' be over 10000 for publication quality results)
#' @return dataframe in which each row gives the statistics (p-value, fold 
#' change and number of standard deviations from the mean) associated with the 
#' enrichment of the stated cell type in the gene list
#' @examples
#' # Load the single cell data
#' ctd <- ctd()
#'
#' # Set the parameters for the analysis
#' # Use 100 bootstrap lists for speed, for publishable analysis use >10000
#' reps <- 100 
#'
#' # Load the data
#' tt_alzh_BA36 <- tt_alzh_BA36()
#' tt_alzh_BA44 <- tt_alzh_BA44()
#'
#' # Run EWCE analysis
#' tt_results_36 <- ewce_expression_data(
#'     sct_data = ctd, tt = tt_alzh_BA36,
#'     annotLevel = 1, ttSpecies = "human", sctSpecies = "mouse"
#' )
#' tt_results_44 <- ewce_expression_data(
#'     sct_data = ctd, tt = tt_alzh_BA44,
#'     annotLevel = 1, ttSpecies = "human", sctSpecies = "mouse"
#' )
#'
#' # Fill a list with the results
#' results <- add.res.to.merging.list(tt_results_36)
#' results <- add.res.to.merging.list(tt_results_44, results)
#'
#' # Perform the merged analysis
#' # For publication reps should be higher
#' merged_res <- merged_ewce(results, reps = 10) 
#' print(merged_res)
#' @export
# @import ggplot2
# @importFrom reshape2 melt
# @importFrom grid unit
# @import plyr
#' @import stats
merged_ewce <- function(results, reps = 100) {
    err_msg <- paste0("ERROR: results list is not valid. Use",
                        " 'add.res.to.merging.list' to generate valid list.")
    # Check arguments are valid
    if (length(results) <= 1) {
        stop(err_msg)
    }

    # If results are directional then seperate directions 
    # and call function recursively
    if (!is.null(results[[1]]$Direction)) {
        for (dir in c("Up", "Down")) {
            found <- 0
            for (i in seq_len(length(results))) {
                if (results[[i]]$Direction == dir) {
                    found <- found + 1
                    if (found == 1) {
                        dir_results <- list(results[[i]])
                    } else {
                        dir_results[[length(dir_results) + 1]] <- results[[i]]
                    }
                    dir_results[[length(dir_results)]]$Direction <- NULL
                }
            }
            if (dir == "Up") {
                merged_res <- 
                    cbind(merged_ewce(dir_results, reps = reps), 
                            Direction = dir)
            } else {
                tmp <- cbind(merged_ewce(dir_results, reps = reps), 
                                Direction = dir)
                merged_res <- rbind(merged_res, tmp)
            }
        }
        return(merged_res)
    } else {
        # If the results are NOT directional
        # First, sum the celltype enrichment values for hit genes in each study
        hit.cells <- results[[1]]$hitCells
        for (i in 2:length(results)) {
            hit.cells <- hit.cells + results[[i]]$hitCells
        }
        # Then:
        # - results[[x]]$bootstrap_data has the enrichment values for the each 
        # of the 10000 original reps
        # - Add these to each other 'reps' times
        count <- 0
        for (i in seq_len(reps)) {
            randomSamples <- 
                sample(seq_len(dim(results[[1]]$bootstrap_data)[1]), 
                        10000, replace = TRUE)
            boot.cells <- results[[1]]$bootstrap_data[randomSamples, ]
            for (i in 2:length(results)) {
                randomSamples <- 
                    sample(seq_len(dim(results[[1]]$bootstrap_data)[1]), 
                            10000, replace = TRUE)
                boot.cells <- 
                    boot.cells + results[[i]]$bootstrap_data[randomSamples, ]
            }
            count <- count + 1
            if (count == 1) {
                all.boot.cells <- boot.cells
            } else {
                all.boot.cells <- rbind(all.boot.cells, boot.cells)
            }
        }
        boot.cells <- all.boot.cells
        p <- fc <- sd_from_mean <- rep(0, length(hit.cells))
        names(p) <- names(fc) <- names(results[[1]]$hitCells)
        for (i in seq_len(length(hit.cells))) {
            p[i] <- 
                sum(hit.cells[i] < boot.cells[, i]) / length(boot.cells[, i])
            fc[i] <- hit.cells[i] / mean(boot.cells[, i])
            sd_from_mean[i] <- 
                (hit.cells[i] - mean(boot.cells[, i])) / sd(boot.cells[, i])
        }
        return(data.frame(CellType = names(p), p = p, fc = fc, 
                            sd_from_mean = sd_from_mean))
    }
}
