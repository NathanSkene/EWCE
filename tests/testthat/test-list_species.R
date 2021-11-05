test_that("list_species works", {
    
    if(!is_32bit()){
        species <- EWCE::list_species()
        testthat::expect_true(methods::is(species,"data.frame"))
        testthat::expect_gte(nrow(species), 21)
    }
})
