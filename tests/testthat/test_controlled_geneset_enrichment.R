# Test for controlled_geneset_enrichment
test_that("Correct controlled & uncontrolled geneset enrichment calculations", {
    # Run test based on vignette data
    data("ctd")
    data("mouse_to_human_homologs")
    m2h <- unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])
    data("schiz_genes")
    mouse.hits.schiz <- unique(m2h[m2h$HGNC.symbol %in% schiz_genes, "MGI.symbol"])
    mouse.bg <- unique(m2h$MGI.symbol)
    data("hpsd_genes")
    mouse.hpsd <- unique(m2h[m2h$HGNC.symbol %in% hpsd_genes, "MGI.symbol"])
    # hyperparameters
    reps <- 100
    # set seed for bootstrap reproducibility - matches vignette
    set.seed(12345678)

    res_hpsd_schiz <- controlled_geneset_enrichment(
        disease_genes = mouse.hits.schiz,
        functional_genes = mouse.hpsd,
        bg_genes = mouse.bg,
        sct_data = ctd, annotLevel = 1,
        reps = reps, controlledCT = "pyramidal CA1"
    )

    # check outputs are as expected
    test1 <- res_hpsd_schiz$reps == reps
    test2 <- res_hpsd_schiz$controlledCT == "pyramidal CA1"
    # from vignette runs we know the following is true - if they change it should throw an erro
    test3 <- res_hpsd_schiz$actualOverlap == 54
    # The probability that functional_genes are enriched in disease_genes
    # while controlling for the level of specificity in controlledCT - p_controlled
    # WITHOUT controlling for the level of specificity in controlledCT - p_uncontrolled
    test4 <- res_hpsd_schiz$p_controlled > res_hpsd_schiz$p_uncontrolled # controlling for cell type worsens prob
    test5 <- res_hpsd_schiz$p_controlled > 0.05 # not sig
    # The z-score that functional_genes are enriched in disease_genes
    # while controlling for the level of specificity in controlledCT - z_controlled
    # WITHOUT controlling for the level of specificity in controlledCT - z_uncontrolled
    test6 <- res_hpsd_schiz$z_controlled < res_hpsd_schiz$z_uncontrolled
    test7 <- res_hpsd_schiz$z_uncontrolled > 1 && res_hpsd_schiz$z_uncontrolled < 2 # z-score between 1 and 2 std above mean
    test8 <- res_hpsd_schiz$z_controlled > 1 && res_hpsd_schiz$z_controlled < 2 # z-score between 1 and 2 std above mean
    # fail if any test doesn't pass
    expect_equal(all(test1, test2, test3, test4, test5, test6, test7, test8), TRUE)
})
