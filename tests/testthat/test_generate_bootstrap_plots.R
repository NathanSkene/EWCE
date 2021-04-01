# test generate_bootstrap_plots and generate_bootstrap_plots_for_transcriptome
test_that("test bootstrap plots and bootstrap plots for transcriptome", {
    # use vignette data
    ctd <- ctd()
    mouse_to_human_homologs <- mouse_to_human_homologs()
    example_genelist <- example_genelist()
    #eh <- query(ExperimentHub::ExperimentHub(), "ewceData")
    #ctd <- eh[["EH5376"]]
    #example_genelist <- eh[["EH5372"]]
    #mouse_to_human_homologs <- eh[["EH5367"]]
    
    m2h <- unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])
    mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% example_genelist, "MGI.symbol"])
    # mouse.bg  = unique(setdiff(m2h$MGI.symbol,mouse.hits))
    mouse.bg <- unique(m2h$MGI.symbol)
    reps <- 10 # <- Use 2 bootstrap lists so it runs quickly, for publishable analysis use >10000
    level <- 1 # <- Use level 1 annotations (i.e. Interneurons)
    full_results <- bootstrap_enrichment_test(
        sct_data = ctd, hits = mouse.hits, bg = mouse.bg,
        reps = reps, annotLevel = level
    )
    options(warn = -1) # turn off warnings for plot warning
    generate_bootstrap_plots(
        sct_data = ctd, hits = mouse.hits,
        bg = mouse.bg, reps = reps, annotLevel = 1,
        full_results = full_results, listFileName = "VignetteGraphs",
        savePath=tempdir()
    )
    options(warn = 0)
    # check the BootstrapPlots folder exists and is non-empty
    #test1 <- dir.exists("~/BootstrapPlots") && length(list.files("~/BootstrapPlots")) > 0
    test1 <- 
        dir.exists(sprintf("%s/BootstrapPlots", tempdir())) && 
        length(list.files(sprintf("%s/BootstrapPlots", tempdir()))) > 0
    # remove folder once tested
    #unlink("~/BootstrapPlots", recursive = TRUE)

    tt_alzh <- tt_alzh()
    #tt_alzh <- eh[["EH5373"]]
    tt_results <- ewce_expression_data(
        sct_data = ctd, tt = tt_alzh, reps=reps,annotLevel = 1,
        ttSpecies = "human", sctSpecies = "mouse"
    )
    options(warn = -1) # turn off warnings for plot warning
    full_results <- generate_bootstrap_plots_for_transcriptome(
        sct_data = ctd, tt = tt_alzh,
        annotLevel = 1, full_results = tt_results,
        listFileName = "examples",
        reps = reps, ttSpecies = "human",
        sctSpecies = "mouse", onlySignif = FALSE
    )
    options(warn = 0)
    # check the BootstrapPlots folder exists and is non-empty
    test2 <- dir.exists("BootstrapPlots") && length(list.files("BootstrapPlots")) > 0
    # remove folder once tested
    unlink("BootstrapPlots", recursive = TRUE)

    # fail if either function didn't create directory and add files
    expect_equal(all(test1, test2), TRUE)
})
