test_that("bootstrap enrichment runs as expected", {
    
    # Test for bootstrap_enrichment_test, using sample data in vignette ensure
    # gives microglia as the only significant enrichment
    if(!is_32bit()){
        set.seed(12345678)
        #### load vignette data ####
        ctd <- ewceData::ctd()
        example_genelist <- ewceData::example_genelist()
        tt_alzh <- ewceData::tt_alzh()
    
        # set input variables
        # Use 10 bootstrap lists so it runs quickly,
        # for publishable analysis use >10000
        reps <- 10
        # Use level 1 annotations (i.e. Interneurons)
        level <- 1
    
        full_results <-
            EWCE::bootstrap_enrichment_test(
                sct_data = ctd,
                sctSpecies = "mouse",
                hits = example_genelist,
                genelistSpecies = "human",
                reps = reps,
                annotLevel = level
            )
        #### Get sig cell types ####
        ewce_sig_cell_types <- as.character(
            full_results$results[full_results$results$p < 0.05, "CellType"]
        )
        # Fail if microglia not returned
        # (occasionally astrocytes will also appear enriched with low reps)
        testthat::expect_true("microglia" %in% ewce_sig_cell_types)
        #### ewce_plot ####
        ewce_plot_res <- EWCE::ewce_plot(
            total_res = full_results$results,
            mtc_method = "BH"
        )$plain
        # Fail if any but ggplot returned
        testthat::expect_true(is(ewce_plot_res, "gg"))
    
        #----------------------------------------------------------
        # Check generate_bootstrap_plots and
        # generate_bootstrap_plots_for_transcriptome
    
        # Use 5 up/down regulated genes (thresh) for speed, default is 250
        thresh <- 5
        options(warn = -1) # turn off warnings for plot warning
        boot_plot_dir1 <- EWCE::generate_bootstrap_plots(
            sct_data = ctd,
            sctSpecies = "mouse",
            hits = example_genelist,
            genelistSpecies = "human",
            annotLevel = level,
            full_results = full_results,
            listFileName = "VignetteGraphs",
            savePath = tempdir(),
            # Important! must match to reps used in bootstrap_enrichment_test
            reps = reps
        )
        options(warn = 0)
        # check the BootstrapPlots folder exists and is non-empty
        testthat::expect_true(
            dir.exists(sprintf("%s/BootstrapPlots", boot_plot_dir1))
        )
        testthat::expect_true(
            length(list.files(sprintf("%s/BootstrapPlots", boot_plot_dir1))) > 0
        )
    
    
        #### tt_results ####
        tt_results <- EWCE::ewce_expression_data(
            sct_data = ctd,
            tt = tt_alzh,
            reps = reps,
            annotLevel = level,
            thresh = thresh,
            ttSpecies = "human",
            sctSpecies = "mouse"
        )
        options(warn = -1) # turn off plot warnings
        boot_plot_dir2 <- EWCE::generate_bootstrap_plots_for_transcriptome(
            sct_data = ctd,
            tt = tt_alzh,
            annotLevel = level,
            full_results = tt_results,
            listFileName = "examples",
            ttSpecies = "human",
            sctSpecies = "mouse",
            onlySignif = FALSE,
            savePath = tempdir(),
            # Important! must match to reps used in ewce_expression_data
            reps = reps
        )
        options(warn = 0)
        # check the BootstrapPlots folder exists and is non-empty
        testthat::expect_true(
            dir.exists(sprintf("%s/BootstrapPlots", boot_plot_dir2))
        )
        testthat::expect_true(
            length(list.files(sprintf("%s/BootstrapPlots", boot_plot_dir2))) > 0
        )
    }
})
