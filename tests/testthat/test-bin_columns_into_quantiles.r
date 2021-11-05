test_that("bin_columns_into_quantiles works", {
    
    if(!is_32bit()){
        set.seed(1234)
        #### Test 1: CTD ####
        ctd <- ewceData::ctd()
        ctd[[1]]$specificity_quantiles <- apply(ctd[[1]]$specificity, 2,
            FUN = bin_columns_into_quantiles,
            numberOfBins = 40
        )
        all_values <- unlist(as.list(ctd[[1]]$specificity_quantiles))
        testthat::expect_equal(sort(unique(all_values)), seq(0, 40))
    
        #### Test 2: When <2 unique non-zero values ####
        mat <- ctd[[1]]$specificity_quantiles
        mat[, 1] <- sample(c(0, 1), size = nrow(mat), replace = TRUE)
        apply
        mat2 <- apply(mat, 2,
            FUN = bin_columns_into_quantiles,
            numberOfBins = 40
        )
        testthat::expect_equal(sort(unique(mat2[, 1])), c(0, 20))
    }
})
