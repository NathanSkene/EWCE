test_that("EWCE expression data creation", {
    # test ewce_expression_data which calls bootstrap_enrichment_test
    # Test for specificity return from generate_celltype_data

    if (!is_32bit()) {
        # bootstrap_enrichment_test tested in depth separately,
        # testing returned specificity here
        set.seed(12345678)
        tt_alzh <- ewceData::tt_alzh()
        ctd <- ewceData::ctd()
        # Use 30 up/down regulated genes (thresh) for speed, default is 250
        thresh <- 30

        tt_results <- EWCE::ewce_expression_data(
            sct_data = ctd,
            tt = tt_alzh,
            reps = 5,
            annotLevel = 1,
            ttSpecies = "human",
            sctSpecies = "mouse",
            thresh = thresh
        )
        up_celltypes <- c(
            names(tt_results$hit.cells.up)[
                tt_results$hit.cells.up == max(tt_results$hit.cells.up)
            ],
            names(tt_results$hit.cells.up)[
                tt_results$hit.cells.up == min(tt_results$hit.cells.up)
            ]
        )
        dwn_celltype <- c(
            names(tt_results$hit.cells.down)[
                tt_results$hit.cells.down == max(tt_results$hit.cells.down)
            ],
            names(tt_results$hit.cells.down)[
                tt_results$hit.cells.down == min(tt_results$hit.cells.down)
            ]
        )

        # check that max and min specificity matches for up and down
        known_maxmin_up_celltypes <- fix_celltype_names(
            c("endothelial-mural", "interneurons")
        )
        known_maxmin_dwn_celltype <- fix_celltype_names(
            c("pyramidal SS", "oligodendrocytes")
        )

        # fail if specificity max and min celltype doesn't match
        testthat::expect_equal(
            (all(up_celltypes == known_maxmin_up_celltypes) &&
                all(dwn_celltype == known_maxmin_dwn_celltype)),
            TRUE
        )

        #----------------------------------------------------------
        # Check add_res_to_merging_list

        # Use 10 up/down regulated genes (thresh) for speed, default is 250
        thresh <- 10
        # For speed set reps to 10, for publication reps should be higher
        reps <- 10

        tt_alzh_BA36 <- ewceData::tt_alzh_BA36()
        tt_alzh_BA44 <- ewceData::tt_alzh_BA44()

        tt_results_36 <- EWCE::ewce_expression_data(
            sct_data = ctd,
            tt = tt_alzh_BA36,
            thresh = thresh,
            reps = reps,
            annotLevel = 1,
            ttSpecies = "human",
            sctSpecies = "mouse"
        )
        tt_results_44 <- EWCE::ewce_expression_data(
            sct_data = ctd,
            tt = tt_alzh_BA44,
            thresh = thresh,
            reps = reps,
            annotLevel = 1,
            ttSpecies = "human",
            sctSpecies = "mouse"
        )
        # Fill a list with the results
        results <- EWCE::add_res_to_merging_list(tt_results)
        results <- EWCE::add_res_to_merging_list(tt_results_36, results)
        results <- EWCE::add_res_to_merging_list(tt_results_44, results)

        # check some non-directional tests
        tt_results_adj <- tt_results[c(1, 2, 4)]
        names(tt_results_adj) <- c(
            "joint_results",
            "hit.cells",
            "bootstrap_data"
        )
        # run the undirectional
        results_adj <- EWCE::add_res_to_merging_list(tt_results_adj)
        # should get error if try to combine directional with undirectional
        error_return <-
            tryCatch(EWCE::add_res_to_merging_list(tt_results_36, results_adj),
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
        test3 <- all.equal(
            results[[1]]$hitCells,
            tt_results$hit.cells.up
        )
        test4 <- all.equal(
            results[[5]]$hitCells,
            tt_results_44$hit.cells.up
        )
        test5 <- all.equal(
            tt_results_36$bootstrap_data.up,
            results[[3]]$bootstrap_data
        )

        # Perform the merged analysis
        merged_res <- EWCE::merged_ewce(results, reps = reps)

        ewce_plot_res <- EWCE::ewce_plot(merged_res)$plain
        # fail if any but ggplot returned
        test6 <- is(ewce_plot_res)[1] == "gg"


        # undirectional tests
        test7 <- length(results_adj[[1]]) == 2
        test8 <- is(error_return, "error")

        # fail if any test doesn't pass
        testthat::expect_true(test1)
        testthat::expect_true(test2)
        testthat::expect_true(test3)
        testthat::expect_true(test4)
        testthat::expect_true(test5)
        testthat::expect_true(test6)
        testthat::expect_true(test7)
        testthat::expect_true(test8)
    }
})
