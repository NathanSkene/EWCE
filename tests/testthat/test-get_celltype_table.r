test_that("get_celltype_table works", {
    if (!is_32bit()) {
        set.seed(1234)
        cortex_mrna <- ewceData::cortex_mrna()

        #### Test 1: get_celltype_table ####
        cortex_mrna$annot$dataset_name <- "cortex_mrna"
        celltype_table <- EWCE::get_celltype_table(cortex_mrna$annot)
        total_celltypes <- length(unique(cortex_mrna$annot$level1class)) +
            length(unique(cortex_mrna$annot$level2class))
        testthat::expect_true(methods::is(celltype_table, "data.frame"))
        testthat::expect_equal(nrow(celltype_table), total_celltypes)
        testthat::expect_equal(ncol(celltype_table), 5)

        #### Test 2: filter_variance_quantiles ####
        exp <- cortex_mrna$exp[seq(1, 300), ]
        ## No normalization
        exp_filt1 <- EWCE:::filter_variance_quantiles(
            exp = exp,
            log10_norm = FALSE
        )
        testthat::expect_equal(nrow(exp_filt1), 180)
        ## SCT normalization
        exp_norm <- EWCE::sct_normalize(exp = exp)
        exp_filt2 <- EWCE:::filter_variance_quantiles(
            exp = exp_norm,
            log10_norm = FALSE
        )
        testthat::expect_equal(nrow(exp_filt2), 180)
        ## Log normalization
        exp_filt3 <- EWCE:::filter_variance_quantiles(
            exp = exp,
            log10_norm = TRUE
        )
        testthat::expect_equal(nrow(exp_filt3), 180)
        ## SCT normalization + Log normalization
        exp_filt4 <- EWCE:::filter_variance_quantiles(
            exp = exp_norm,
            log10_norm = TRUE
        )
        testthat::expect_equal(nrow(exp_filt4), 180)

        ## CONCLUSION:
        ## log normalisation has no effect using computing quantiles 
        ## using the "EWCE" method (default), 
        ## but is essential when computing quantiles using 
        ## the stats::ecdf method.
    }
})
