test_that("check_percent_hits works", {
    if (!is_32bit()) {
        full_results <- EWCE::example_bootstrap_results()
        testthat::expect_true(length(full_results) %in% c(3,4) )
        testthat::expect_true(methods::is(full_results$results, "data.frame"))

        report <- EWCE::check_percent_hits(
            full_results = full_results,
            target_celltype = "microglia"
        )
        testthat::expect_true(
            all(c("target_hits", "percent_hits", "target_celltype")
            %in% names(report))
        )
        testthat::expect_equal(report$percent_hits, 14.3)
    }
})
