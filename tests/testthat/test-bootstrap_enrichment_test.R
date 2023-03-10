test_that("bootstrap enrichment runs as expected", {

    # Test for bootstrap_enrichment_test, using sample data in vignette ensure
    # gives microglia as the only significant enrichment
    if (!is_32bit()) {
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

        #----------------------------------------------------------
        # Check generate_bootstrap_plots and
        # generate_bootstrap_plots_for_transcriptome

        # Use 5 up/down regulated genes (thresh) for speed, default is 250
        thresh <- 5
        options(warn = -1) # turn off warnings for plot warning
        boot_out1 <- generate_bootstrap_plots(
            sct_data = ctd,
            sctSpecies = "mouse",
            hits = example_genelist,
            genelistSpecies = "human",
            annotLevel = level,
            full_results = full_results,
            listFileName = "VignetteGraphs",
            # Important! must match to reps used in bootstrap_enrichment_test
            reps = reps
        )
        options(warn = 0) 
        testthat::expect_length(boot_out1$paths,4)
        for(f in boot_out1$paths){
            testthat::expect_true(
                file.exists(f)
            ) 
        }
        testthat::expect_length(boot_out1$plots,4)
        for(p in boot_out1$plots){
            testthat::expect_true(
                methods::is(p,"gg")
            ) 
        }  
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
        boot_dir2 <- generate_bootstrap_plots_for_transcriptome(
            sct_data = ctd,
            tt = tt_alzh,
            annotLevel = level,
            full_results = tt_results,
            listFileName = "examples",
            ttSpecies = "human",
            sctSpecies = "mouse",
            sig_only = FALSE, 
            # Important! must match to reps used in ewce_expression_data
            reps = reps,
            # Only do one plot type for demo purposes
            # (other plot types take absurdly long)
            plot_types = "bootstrap"
        )
        options(warn = 0)
        # check the BootstrapPlots folder exists and is non-empty 
        testthat::expect_true(
            dir.exists(boot_dir2)
        )
        testthat::expect_true(
            length(list.files(boot_dir2)) > 0
        )
    }
})
