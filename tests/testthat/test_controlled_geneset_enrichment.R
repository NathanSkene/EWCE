test_that("Correct controlled & uncontrolled geneset enrichment calculations", {
    # Test for controlled_geneset_enrichment
    
    if(!is_32bit()){
        # Run test based on vignette data
        ctd <- ewceData::ctd()
        mouse_to_human_homologs <- ewceData::mouse_to_human_homologs()
    
        m2h <- unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])
        schiz_genes <- ewceData::schiz_genes()
        mouse.hits.schiz <- unique(m2h[m2h$HGNC.symbol %in% schiz_genes, "MGI.symbol"])
        mouse.bg <- unique(m2h$MGI.symbol)
        hpsd_genes <- ewceData::hpsd_genes()
        mouse.hpsd <- unique(m2h[m2h$HGNC.symbol %in% hpsd_genes, "MGI.symbol"])
        # hyperparameters
        reps <- 100
        # set seed for bootstrap reproducibility - matches vignette
        set.seed(12345678)
    
        res_hpsd_schiz <- controlled_geneset_enrichment(
            disease_genes = mouse.hits.schiz,
            functional_genes = mouse.hpsd,
            disease_genes_species = "mouse",
            functional_genes_species = "mouse",
            sct_data = ctd,
            # bg = mouse.bg,
            sctSpecies = "mouse",
            output_species = "mouse",
            annotLevel = 1,
            reps = reps,
            controlledCT = "pyramidal CA1"
        )
        ### Input as original human gene lists ####
        # res_hpsd_schiz2 <- EWCE::controlled_geneset_enrichment(
        #     disease_genes = schiz_genes,
        #     functional_genes = hpsd_genes,
        #     disease_genes_species = "human",
        #     functional_genes_species = "human",
        #     sct_data = ctd,
        #     sctSpecies = "mouse",
        #     output_species = "human",
        #     annotLevel = 1,
        #     reps = reps,
        #     controlledCT = "pyramidal CA1"
        # )
    
        #### Check args are as expected ####
        testthat::expect_true(res_hpsd_schiz$reps == reps)
        testthat::expect_true(res_hpsd_schiz$controlledCT == "pyramidal CA1")
    
        #### from vignette runs we know the following is true ####
        # if they change it should throw an error
        testthat::expect_true(res_hpsd_schiz$actualOverlap == 54)
    
        # The probability that functional_genes are enriched in disease_genes
        # while controlling for the level of specificity
        # in controlledCT - p_controlled
        # WITHOUT controlling for the level of specificity
        # in controlledCT - p_uncontrolled
        # controlling for cell type worsens prob
        testthat::expect_gte(
            res_hpsd_schiz$p_controlled,
            res_hpsd_schiz$p_uncontrolled
        )
    
        #### Check that the results are not significant ####
        ## May vary depending on background used ####
        # testthat::expect_true(res_hpsd_schiz$p_controlled > 0.05)
    
        # The z-score that functional_genes are enriched in disease_genes
        # while controlling for the level of specificity in controlledCT
        #- z_controlled
        # WITHOUT controlling for the level of specificity in controlledCT
        #- z_uncontrolled
        testthat::expect_lte(
            res_hpsd_schiz$z_controlled,
            res_hpsd_schiz$z_uncontrolled
        )
    
        #### z-score between 1 and 2 std above mean ####
        # testthat::expect_true(res_hpsd_schiz$z_uncontrolled > 1 &&
        #                           res_hpsd_schiz$z_uncontrolled < 2)
    
        #### z-score between 1 and 2 std above mean ####
        testthat::expect_true(res_hpsd_schiz$z_controlled > 1 &&
            res_hpsd_schiz$z_controlled < 2)
    }
})
