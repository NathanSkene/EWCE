#' Add to results to merging list
#'
#' \code{add_res_to_merging_list} Adds EWCE results to a list
#' for merging analysis
#'
#' @param full_res results list generated using
#' \code{\link{bootstrap_enrichment_test}} or \code{\link{ewce_expression_data}}
#' functions. Multiple results tables can be merged into one
#' results table, as long as the 'list' column is set to distinguish them.
#' @param existing_results Output of previous rounds from adding results to
#' list. Leave empty if this is the first item in the list.
#' @return merged results list
#' @examples
#' # Load the single cell data
#' library(ewceData)
#' ctd <- ctd()
#'
#' # Load the data
#' tt_alzh <- tt_alzh()
#' #tt_alzh_BA36 <- tt_alzh_BA36()
#' # Use 3 bootstrap lists for speed, for publishable analysis use >10000
#' reps <- 3
#' # Use 5 up/down regulated genes (thresh) for speed, default is 250
#' thresh = 5
#' # Run EWCE analysis
#' #tt_results <- ewce_expression_data(
#' #    sct_data = ctd, tt = tt_alzh, annotLevel = 1, thresh = thresh,
#' #    reps = reps, ttSpecies = "human", sctSpecies = "mouse"
#' #)
#' #tt_results_36 <- ewce_expression_data(
#' #    sct_data = ctd, tt = tt_alzh_BA36, annotLevel = 1, thresh = thresh,
#' #    reps = reps, ttSpecies = "human", sctSpecies = "mouse"
#' #)
#'
#' # Fill a list with the results
#' results <- add_res_to_merging_list(tt_alzh)
#' #results <- add_res_to_merging_list(tt_alzh_BA36, results)
#' @export
add_res_to_merging_list <- function(full_res, existing_results = NULL) {
    msg <- "ERROR: Cannot merge directional with non-directional results table"
    # Check if the results set is directional
    if (is.null(full_res$bootstrap_data)) {
        # If it is directional, and other results are in list,
        # check they are all directional
        if (!is.null(existing_results)) {
            for (i in seq_along(length(existing_results))) {
                if (is.null(existing_results[[i]]$Direction)) {
                    stop(msg)
                }
            }
        }
        # Setup the new entries
        new_entry_up <- list(
            Direction = "Up", bootstrap_data = full_res$bootstrap_data.up,
            hitCells = full_res$hit.cells.up
        )
        new_entry_down <- list(
            Direction = "Down", bootstrap_data = full_res$bootstrap_data.down,
            hitCells = full_res$hit.cells.down
        )
        # Add them to the list
        if (is.null(existing_results)) {
            existing_results <- list(new_entry_up)
        } else {
            existing_results[[length(existing_results) + 1]] <- new_entry_up
        }
        existing_results[[length(existing_results) + 1]] <- new_entry_down
    } else {
        # If the list is NOT directional
        # If it is NOT directional, and other results are in list,
        # check they are all not directional
        if (!is.null(existing_results)) {
            for (i in seq_along(length(existing_results))) {
                if (!is.null(existing_results[[i]]$Direction)) {
                    stop(msg)
                }
            }
        }
        # Setup the new entries
        new_entry <- list(
            bootstrap_data = full_res$bootstrap_data,
            hitCells = full_res$hit.cells
        )
        # Add new entry to the list
        if (is.null(existing_results)) {
            existing_results <- list(new_entry)
        } else {
            existing_results[[length(existing_results) + 1]] <- new_entry
        }
    }
    return(existing_results)
}
