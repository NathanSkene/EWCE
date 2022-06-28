test_that("merge_sce works", {
    if (!is_32bit()) {
        ctd <- ewceData::ctd()

        #### Test 1: convert CTD to SCE list ####
        sce_list <- EWCE::ctd_to_sce(object = ctd)
        testthat::expect_true(methods::is(sce_list, "list"))
        testthat::expect_length(sce_list, length(ctd))
        testthat::expect_true(methods::is(sce_list$level_1, "SingleCellExperiment"))

        ### Test 2: merge sce_list ####
        sce_combine <- EWCE::merge_sce(sce_list = sce_list)
        testthat::expect_true(methods::is(sce_combine, "SingleCellExperiment"))
        testthat::expect_gte(nrow(sce_combine), 15000)
        testthat::expect_equal(ncol(sce_combine), 55)

        #### Test 3: merge CTD ####
        ## Let's pretend these are different CTD datasets
        ctd2 <- ctd1 <- ctd
        CTD_list <- list(ctd1, ctd2)
        SCE_merged <- EWCE::merge_ctd(
            CTD_list = CTD_list,
            as_SCE = TRUE,
            gene_union = TRUE
        )
        testthat::expect_length(SCE_merged, length(CTD_list[[1]]))
        testthat::expect_false(EWCE:::is_celltypedataset(SCE_merged$level_1))
        testthat::expect_true(
            methods::is(SCE_merged$level_1, "SingleCellExperiment")
        )


        #### Test 4: convert between old/new CTD formats ####
        ## New to old
        # Level 1
        ctd_old_lvl1 <- EWCE:::convert_new_ewce_to_old(
            ctd = ctd,
            lvl = 1
        )
        testthat::expect_true(
            all(c("cell_dists", "all_scts") %in% names(ctd_old_lvl1[[1]]))
        )
        lvl1_path <- tempfile(fileext = ".rda")
        save(ctd_old_lvl1, file = lvl1_path)
        # Level 2
        ctd_old_lvl2 <- EWCE:::convert_new_ewce_to_old(
            ctd = ctd,
            lvl = 2
        )
        testthat::expect_true(
            all(c("cell_dists", "all_scts") %in% names(ctd_old_lvl2[[1]]))
        )
        lvl2_path <- tempfile(fileext = ".rda")
        save(ctd_old_lvl2, file = lvl2_path)
        ## Old to new
        # Method 1
        ctd_new1 <- EWCE:::convert_old_ewce_to_new(
            level1 = lvl1_path,
            level2 = lvl2_path
        )
        testthat::expect_length(ctd_new1, 2)
        testthat::expect_true(EWCE:::is_celltypedataset(ctd_new1))
        # Method 2
        ctd_new2 <- EWCE:::convert_old_ewce_to_new(
            celltype_data = list(
                ctd_old_lvl1[[1]],
                ctd_old_lvl2[[1]]
            )
        )
        testthat::expect_length(ctd_new2, 2)
        testthat::expect_true(EWCE:::is_celltypedataset(ctd_new2))

        #### Test 5: plot_ctd ####
        gp <- EWCE::plot_ctd(
            ctd = ctd,
            genes = sample(rownames(ctd[[1]]$mean_exp), 4)
        )
        testthat::expect_true(methods::is(gp, "gg"))
        
        
        #### Test 6: merge CTD_list into one CTD 
        CTD_merged <- EWCE::merge_ctd(CTD_list = CTD_list)
        testthat::expect_true(EWCE:::is_celltypedataset(CTD_merged))
        testthat::expect_equal(length(CTD_merged), length(CTD_list[[1]]))
    }
})
