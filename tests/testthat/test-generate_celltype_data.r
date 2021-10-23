test_that("generate_celltype_data works", {

    # Load the single cell data
    cortex_mrna <- ewceData::cortex_mrna()
    # Use only a subset to keep the example quick
    expData <- cortex_mrna$exp[seq(1, 100), ]
    l1 <- cortex_mrna$annot$level1class
    l2 <- cortex_mrna$annot$level2class
    annotLevels <- list(l1 = l1, l2 = l2)
    #### As DelayedArray ####
    res <- EWCE::generate_celltype_data(
        exp = expData,
        annotLevels = annotLevels,
        convert_orths = TRUE,
        input_species = "mouse",
        output_species = "human",
        # Converts expData to DelayedArray before processed further.
        # Doesn't convert CTD matrices into DelayedArray.
        as_DelayedArray = TRUE,
        groupName = "allKImouse",
        return_ctd = TRUE
    )
    testthat::expect_true(EWCE:::is_celltypedataset(res$ctd))
    testthat::expect_true(EWCE:::is_sparse_matrix(res$ctd[[1]]$mean_exp))
    testthat::expect_true(file.exists(res$fNames))
})
