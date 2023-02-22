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
})
