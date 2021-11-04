test_that("list_species works", {
    species <- EWCE::list_species()
    testthat::expect_true(methods::is(species,"data.frame"))
    testthat::expect_gte(nrow(species), 21)
})
