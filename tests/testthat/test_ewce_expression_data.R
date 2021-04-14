# test ewce_expression_data which calls bootstrap_enrichment_test
# Test for specificity return from generate_celltype_data
test_that("EWCE expression data creation", {
    # bootstrap_enrichment_test tested in depth separately, testing returned specificity here
    set.seed(12345678)
    tt_alzh <- tt_alzh()
    ctd <- ctd()
    #eh <- query(ExperimentHub::ExperimentHub(), "ewceData")
    #tt_alzh <- eh[["EH5373"]]
    #ctd <- eh[["EH5376"]]

    # Use 30 up/down regulated genes (thresh) for speed, default is 250
    thresh = 30

    tt_results <- ewce_expression_data(sct_data = ctd, tt = tt_alzh, reps = 5,
                                        annotLevel = 1, ttSpecies = "human",
                                        sctSpecies = "mouse", thresh=thresh)
    up_celltypes <- c(
        names(tt_results$hit.cells.up)[tt_results$hit.cells.up == max(tt_results$hit.cells.up)],
        names(tt_results$hit.cells.up)[tt_results$hit.cells.up == min(tt_results$hit.cells.up)]
    )
    dwn_celltypes <- c(
        names(tt_results$hit.cells.down)[tt_results$hit.cells.down == max(tt_results$hit.cells.down)],
        names(tt_results$hit.cells.down)[tt_results$hit.cells.down == min(tt_results$hit.cells.down)]
    )

    # check that max and min specificity matches for up and down
    known_maxmin_up_celltypes <- c("endothelial-mural", "interneurons")
    known_maxmin_dwn_celltypes <- c("pyramidal SS", "oligodendrocytes")

    # fail if specificity max and min cell types doesn't match
    expect_equal(
        (all(up_celltypes == known_maxmin_up_celltypes) &&
            all(dwn_celltypes == known_maxmin_dwn_celltypes)),
        TRUE
    )
    
    #----------------------------------------------------------
    #Check add_res_to_merging_list
    
    # Use 10 up/down regulated genes (thresh) for speed, default is 250
    thresh = 10
    # For speed set reps to 10, for publication reps should be higher
    reps = 10
    
    tt_alzh_BA36 <- tt_alzh_BA36()
    tt_alzh_BA44 <- tt_alzh_BA44()
    
    tt_results_36 <- ewce_expression_data(
        sct_data = ctd, tt = tt_alzh_BA36,thresh=thresh, reps=reps,
        annotLevel = 1, ttSpecies = "human",
        sctSpecies = "mouse"
    )
    tt_results_44 <- ewce_expression_data(
        sct_data = ctd, tt = tt_alzh_BA44,thresh=thresh, reps=reps,
        annotLevel = 1, ttSpecies = "human",
        sctSpecies = "mouse"
    )
    # Fill a list with the results
    results <- add_res_to_merging_list(tt_results)
    results <- add_res_to_merging_list(tt_results_36, results)
    results <- add_res_to_merging_list(tt_results_44, results)
    
    # check some non-directional tests
    tt_results_adj <- tt_results[c(1, 2, 4)]
    names(tt_results_adj) <- c("joint_results", "hit.cells", "bootstrap_data")
    # run the undirectional
    results_adj <- add_res_to_merging_list(tt_results_adj)
    # should get error if try to combine directional with undirectional
    error_return <-
        tryCatch(add_res_to_merging_list(tt_results_36, results_adj),
                 error = function(e) e,
                 warning = function(w) w
        )
    
    # check returns
    test1 <- length(results) == 6
    # make sure all metrics returned for each
    test2 <- all.equal(
        unlist(lapply(results, function(x) names(x))),
        rep(c("Direction", "bootstrap_data", "hitCells"), 6)
    )
    # check output equal to input before combination
    test3 <- all.equal(results[[1]]$hitCells, tt_results$hit.cells.up)
    test4 <- all.equal(results[[5]]$hitCells, tt_results_44$hit.cells.up)
    test5 <- all.equal(tt_results_36$bootstrap_data.up, results[[3]]$bootstrap_data)
    
    # Perform the merged analysis
    merged_res <- merged_ewce(results, reps = reps)
    
    ewce_plot_res <- ewce_plot(merged_res)$plain
    # fail if any but ggplot returned
    test6 <- is(ewce_plot_res)[1] == "gg"
    
    
    # undirectional tests
    test7 <- length(results_adj[[1]]) == 2
    test8 <- is(error_return, "error")
    
    # fail if any subtest isn't true
    expect_equal(all(test1, test2, test3, test4, test5, test6, test7, test8), TRUE)
})
