test_that("Filter genes without 1 to 1 homolog test", {

    # Test for filter_genes_without_1to1_homolog, using sample data in vignette
    # inspect drop_uninformative_genes
    # Use Vignette Dataset to check function output
    if (!is_32bit()) {
        cortex_mrna <- ewceData::cortex_mrna()
        if (!file.exists(sprintf("%s/MRK_List2.rpt", tempdir()))) {
            download.file(
                "http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt",
                destfile = sprintf("%s/MRK_List2.rpt", tempdir())
            )
        }
        nKeep <- 1000
        must_keep <- c("Apoe", "Gfap", "Gapdh")
        set.seed(123458)
        keep_genes <- c(must_keep, sample(rownames(cortex_mrna$exp), 997))
        cortex_mrna$exp <- cortex_mrna$exp[keep_genes, ]

        # NOTE not normalising for this test, too time intensive
        exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(
            exp = cortex_mrna$exp,
            ## Prevent from converting genes to human
            ## (thus dropping more than expected)
            input_species = "mouse",
            output_species = "mouse",
            level2annot = cortex_mrna$annot$level2class,
            drop_nonhuman_genes = FALSE,
            ### IMPORTANT! Testthat causes conflicts with parallelization
            no_cores = 1
        )
        # check the number of dropped genes is the same
        # as expected from past runs: 248 else fail
        testthat::expect_equal(
            nrow(cortex_mrna$exp) - nrow(exp_CortexOnly_DROPPED),
            248
        )
        #----------------------------------------------------------
        #### Check filter_genes_without_1to1_homolog ###
        annotLevels <- list(
            level1class = cortex_mrna$annot$level1class,
            level2class = cortex_mrna$annot$level2class
        )
        fNames_CortexOnly <- generate_celltype_data(
            exp = exp_CortexOnly_DROPPED,
            annotLevels = annotLevels,
            groupName = "kiCortexOnly",
            input_species = "mouse",
            output_species = "mouse",
            savePath = tempdir(),
            ### IMPORTANT! Testthat causes conflicts with parallelization
            no_cores = 1
        )
        #### filter  orthologs ####
        # Deprecated function
        testthat::expect_warning(
            fnames_warn <- EWCE:::filter_genes_without_1to1_homolog(
                filenames = fNames_CortexOnly
            )
        )
        # Current function
        fNames_CortexOnly <- EWCE::filter_nonorthologs(
            filenames = fNames_CortexOnly
        )
        # load and inspect that none were removed
        ctd <- EWCE::load_rdata(fNames_CortexOnly[1])
        # Check that genes were dropped
        testthat::expect_lte(
            nrow(ctd[[1]]$specificity),
            nrow(exp_CortexOnly_DROPPED)
        )
        # check the number of dropped genes is
        # the same as expected from past runs
        testthat::expect_gte(
            nrow(ctd[[1]]$specificity),
            650
        )
    }
})
