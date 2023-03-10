#' check_args_for_bootstrap_plot_generation
#' 
#' Check the input arguments of the 
#' \link[EWCE]{generate_bootstrap_plots_for_transcriptome}.
#' 
#' @inheritParams generate_bootstrap_plots_for_transcriptome
#' @return Null output.
#' 
#' @keywords internal
check_args_for_bootstrap_plot_generation <- function(sct_data, 
                                                     tt, 
                                                     thresh,
                                                     annotLevel, 
                                                     reps,
                                                     full_results,
                                                     listFileName,
                                                     showGNameThresh, 
                                                     sortBy) {
    # Check the arguments
    if(all(is.na(full_results))) {
        stop_msg <- paste("Must provide valid full_results",
                          "from ewce_expression_data function.")
        stop(stop_msg)
    }
    correct_length <- length(full_results) == 5
    required_names <- c(
        "joint_results", "hit.cells.up",
        "hit.cells.down", "bootstrap_data.up",
        "bootstrap_data.down"
    )
    all_required_names <- sum(names(full_results) %in% required_names) == 5
    err_msg <- paste0(
        "ERROR: full_results is not valid output from the",
        " ewce_expression_data function. This function only",
        " takes data generated from transcriptome analyses."
    )
    if (!correct_length | !all_required_names) {
        stop(err_msg)
    }

    # Check the arguments
    err_msg2 <- paste0(
        "ERROR: tt does not contain a column with value",
        " passed in sortBy argument"
    )
    if (!sortBy %in% colnames(tt)) {
        stop(err_msg2)
    }
    err_msg3 <- paste0(
        "ERROR: length of table is less than twice the",
        " size of threshold"
    )
    if (dim(tt)[1] < (thresh * 2)) {
        stop(err_msg3)
    } 
}
