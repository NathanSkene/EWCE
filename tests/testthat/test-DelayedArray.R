test_that("DelayedArray works", {

    if(!is_32bit()){
        #### Setup data ####
        cortex_mrna <- ewceData::cortex_mrna()
        expMatrix <- DelayedArray::DelayedArray(cortex_mrna$exp)
        # ctd <- ewceData::ctd()
        ctd <- list(level1 = list(),
                    level2 = list())
        ctd[[1]][["annot"]] <- cortex_mrna$annot$level1class
        #### Set DelayedArray parameters ####
        EWCE:::assign_cores(worker_cores = 1)
        DelayedArray:::set_verbose_block_processing(verbose = TRUE)
    
        #### Test calculate_meanexp_for_level ####
        ctd_oneLevel <- calculate_meanexp_for_level(
            ctd_oneLevel = ctd[[1]],
            expMatrix = expMatrix
        )
        testthat::expect_length(ctd_oneLevel, 2)
        testthat::expect_true(all(c("annot", "mean_exp") %in% names(ctd_oneLevel)))
    
        #### Test calculate_specificity_for_level ####
        ctd_oneLevel_mod <- calculate_specificity_for_level(
            ctd_oneLevel = ctd_oneLevel
        )
        testthat::expect_length(ctd_oneLevel_mod, 3)
        testthat::expect_true(
            all(c("annot", "mean_exp", "specificity") %in% names(ctd_oneLevel_mod))
        )
        
        #### Test delayedarray_normalize ####
        exp_norm <- EWCE:::delayedarray_normalize(exp = expMatrix,
                                                  log_norm = TRUE, 
                                                  min_max = TRUE)
        testthat::expect_equal(dim(exp_norm), dim(expMatrix)) 
    }
})
