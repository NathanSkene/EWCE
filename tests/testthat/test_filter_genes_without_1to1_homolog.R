# Test for filter_genes_without_1to1_homolog, using sample data in vignette
test_that("Filter genes without 1 to 1 homolog test", {
    # inspect drop_uninformative_genes
    # Use Vignette Dataset to check function output
    cortex_mrna <- cortex_mrna()
    #eh <- query(ExperimentHub::ExperimentHub(), "ewceData")
    #cortex_mrna <- eh[["EH5381"]]
    
    #if (!file.exists("MRK_List2.rpt")) {
    if (!file.exists(sprintf("%s/MRK_List2.rpt", tempdir()))){
        download.file("http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt", 
                        destfile = sprintf("%s/MRK_List2.rpt", tempdir()))
    }
    nKeep <- 1000
    must_keep <- c("Apoe", "Gfap", "Gapdh")
    set.seed(123458)
    keep_genes <- c(must_keep, sample(rownames(cortex_mrna$exp), 997))
    cortex_mrna$exp <- cortex_mrna$exp[keep_genes, ]

    # NOTE not normalising for this test, too time intensive

    exp_CortexOnly_DROPPED <- drop_uninformative_genes(
        exp = cortex_mrna$exp,
        level2annot = cortex_mrna$annot$level2class
    )
    
    # check the number of dropped genes is the same as expected from past runs: 248 else fail
    expect_equal(nrow(cortex_mrna$exp) - nrow(exp_CortexOnly_DROPPED), 248)
    
    #----------------------------------------------------------
    #Check filter_genes_without_1to1_homolog
    
    annotLevels <- list(level1class = cortex_mrna$annot$level1class, level2class = cortex_mrna$annot$level2class)
    fNames_CortexOnly <- generate_celltype_data(exp = exp_CortexOnly_DROPPED, 
                                                    annotLevels = annotLevels, 
                                                    groupName = "kiCortexOnly",
                                                    savePath = tempdir())

    # filter only orthologs
    fNames_CortexOnly <- filter_genes_without_1to1_homolog(fNames_CortexOnly)
    # load and inspect that none were removed
    load(fNames_CortexOnly[1])
    nrow(ctd[[1]]$specificity)

    # remove folder once tested
   # unlink("MRK_List2.rpt", recursive = TRUE)

    # check the number of dropped genes is the same as expected from past runs: 248 else fail
    expect_equal(nrow(ctd[[1]]$specificity), nrow(exp_CortexOnly_DROPPED))
})
