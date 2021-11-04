test_that("merge_sce works", {
    set.seed(1234)
    ctd <- ewceData::ctd()
    
    #### Test 1: convert CTD to SCE list #### 
    sce_list <- EWCE::ctd_to_sce(object = ctd)
    testthat::expect_true(methods::is(sce_list,"list"))
    testthat::expect_length(sce_list, length(ctd))
    testthat::expect_true(methods::is(sce_list$level_1,"SingleCellExperiment"))
    
    ### Test 2: merge sce_list ####
    sce_combine <- EWCE::merge_sce(sce_list = sce_list) 
    total_celltypes <- length(unique(ctd[[1]]$annot)) +  
                       length(unique(ctd[[2]]$annot))
    testthat::expect_true(methods::is(sce_combine,"SingleCellExperiment"))
    testthat::expect_equal(nrow(sce_combine), nrow(ctd[[1]]$mean_exp))
    testthat::expect_equal(ncol(sce_combine), total_celltypes)
    
    #### Test 3: merge CTD ####
    ## Let's pretend these are different CTD datasets 
    ctd2 <- ctd1 <- ctd
    CTD_list <- list(ctd1, ctd2)
    SCE_merged <- EWCE::merge_ctd(
        CTD_list = CTD_list,
        as_SCE = TRUE,
        gene_union = TRUE
    )
    testthat::expect_length(SCE_merged, length(CTD_list))
    testthat::expect_length(SCE_merged, length(CTD_list[[1]]))
    testthat::expect_true(methods::is(SCE_merged$level_1,
                                      "SingleCellExperiment"))
    
    #### Test 4: is_celltypedataset ####
    testthat::expect_true(EWCE:::is_celltypedataset(ctd))
    testthat::expect_false(EWCE:::is_celltypedataset(CTD_list))
    
    #### Test 5: plot_ctd ####
    gp <- EWCE::plot_ctd(ctd = ctd,
                         genes = sample(rownames(ctd[[1]]$mean_exp),4))
    testthat::expect_true(methods::is(gp,"gg"))
})
