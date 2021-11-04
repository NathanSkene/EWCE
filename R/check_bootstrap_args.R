check_bootstrap_args <- function(sct_data,
    hits,
    annotLevel,
    reps,
    controlledCT = NULL) {
    #### Check an SCT dataset was provided ####
    if (unique(is.null(sct_data)) ||
        (!is_celltypedataset(ctd = sct_data))) {
        stop("Must provide valid CellTypeDataset to sct_data.")
    }
    #### Check an hits was provided ####
    if ((!exists("hits")) || unique(is.null(hits))) {
        stop("Must provide a gene list to hits.")
    }
    #### Check an hits was provided ####
    if ((!exists("annotLevel")) ||
        (annotLevel < 1) ||
        (annotLevel > length(sct_data))) {
        stop("Must provide a valid annotLevel <= ", length(sct_data), ".")
    }
    if ((!exists("reps")) ||
        (reps < 1)) {
        stop("Must provide a valid reps > 0 .")
    }
    #### Check if controlling for another celltype ###
    if (!is.null(controlledCT)) {
        if (!controlledCT %in% colnames(sct_data[[1]]$specificity)) {
            err_msg <- paste0(
                "invalid celltype name passed in controlledCT.",
                " This argument is optional. Leave empty if you do not",
                " wish to control for a celltypes expression."
            )
            stop(err_msg)
        }
    }
}
