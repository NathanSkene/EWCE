test_that("filter_ctd_genes works", {
 
    ctd <- ewceData::ctd() 
    n_genes <- 100
    gene_subset <- rownames(ctd[[1]]$mean_exp)[seq_len(n_genes)]
    
    #### Works with old CTD format ####
    ctd_subset <- filter_ctd_genes(ctd = ctd,
                                   gene_subset = gene_subset)
    testthat::expect_true(
        all(lapply(ctd_subset, function(x)nrow(x$mean_exp))==n_genes)
    )
    
    #### Works with new CTD format ####
    ctd <- standardise_ctd(ctd, input_species="mouse")
    gene_subset <- rownames(ctd[[1]]$mean_exp)[seq_len(n_genes)]
    ctd_subset <- filter_ctd_genes(ctd = ctd,
                                   gene_subset = gene_subset)
    testthat::expect_true(
        all(lapply(ctd_subset, function(x)nrow(x$mean_exp))==n_genes)
    )
})
