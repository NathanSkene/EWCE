check_args_for_bootstrap_plot_generation <- function(sct_data, tt, thresh,
    annotLevel, reps,
    full_results,
    listFileName,
    showGNameThresh,
    ttSpecies, sctSpecies,
    sortBy) {
    # Check the arguments
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

    # Check that the top table has correct columns
    err_msg4 <- paste0(
        "ERROR: if ttSpecies==human then there must be an",
        " HGNC.symbol column"
    )
    err_msg5 <- paste0(
        "ERROR: if ttSpecies==human then there must be an",
        " MGI.symbol column"
    )
    if (ttSpecies == "human" & !"HGNC.symbol" %in% colnames(tt)) {
        stop(err_msg4)
    }
    if (ttSpecies == "mouse" & !"MGI.symbol" %in% colnames(tt)) {
        stop(err_msg5)
    }

    if (ttSpecies == "human" & sctSpecies == "human") {
        tt$MGI.symbol <- tt$HGNC.symbol
    }
    mouse_to_human_homologs <- ewceData::mouse_to_human_homologs()
    m2h <- mouse_to_human_homologs[, c("MGI.symbol", "HGNC.symbol")]
    if (ttSpecies == "human" & sctSpecies == "mouse") {
        tt <- merge(tt, m2h, by = "HGNC.symbol")
    }
    if (ttSpecies == "mouse" & sctSpecies == "human") {
        tt <- merge(tt, m2h, by = "HGNC.symbol")
        tt$MGI.symbol <- tt$HGNC.symbol
    }
    return(tt)
}
