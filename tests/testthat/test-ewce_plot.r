test_that("ewce_plot works", {
 
    full_results <- EWCE::example_bootstrap_results()
    ctd <- ewceData::ctd()
    #### ewce_plot ####
    ewce_plot_res <- ewce_plot(
        total_res = full_results$results, 
        ctd = ctd, 
        make_dendro = TRUE
    )
    # Fail if any but ggplot returned
    testthat::expect_true(methods::is(ewce_plot_res$plain, "gg"))
    testthat::expect_true(methods::is(ewce_plot_res$withDendro, "gg"))
    #remove short cut dendrogram res from ctd anbd rerun, should get same order
    ctd_basic <- ctd
    ctd_basic[[1]]$plotting <- NULL
    #### ewce_plot ####
    ewce_plot_res_basic <- ewce_plot(
      total_res = full_results$results, 
      ctd = ctd_basic, 
      make_dendro = TRUE
    )
    #so order of 4 plots should be the same
    testthat::expect_true(
      all(sapply(list(ewce_plot_res_basic$withDendro$data$CellType, 
                      ewce_plot_res$plain$data$CellType, 
                      ewce_plot_res_basic$plain$data$CellType), 
                 FUN = identical, ewce_plot_res$withDendro$data$CellType))
    )
})
