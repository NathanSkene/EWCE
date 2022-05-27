test_that("method to remove/fix an expected set of genes", {
    # Test for fix_bad_mgi_symbols

    if (!is_32bit()) {
        requireNamespace("data.table")
        # Use Vignette Dataset to check function, alter input gene names
        cortex_mrna <- ewceData::cortex_mrna()

        # IMPORTANT!: Do not delete. Check that outputs "EWCE_return" will fail otherwise.
        if (!file.exists(sprintf("%s/MRK_List2.rpt", tempdir()))) {
            options(timeout = 60*3) # huge file, takes a long time
            utils::download.file(
                "http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt",
                destfile = sprintf("%s/MRK_List2.rpt", tempdir())
            )
        }
        # take a subset for testing
        # causes error when no mismatch
        test_exp_set <- cortex_mrna$exp[seq(1, 10), seq(1, 50)]
        # Add in fake gene data for which the gene could have issues with excel:
        # MARCH1 gene would go to Mar-01
        rownames(test_exp_set)[7] <- "Mar-01"

        # catch error when no exp inputted
        error_return <-
            tryCatch(EWCE::fix_bad_mgi_symbols(
                exp = NULL # ,
                # mrk_file_path = sprintf("%s/MRK_List2.rpt", tempdir())
            ),
            error = function(e) e,
            warning = function(w) w
            )

        # catch error when no test_exp_set inputted as character
        test_exp_set_char <- apply(test_exp_set, 2, as.character)
        rownames(test_exp_set_char) <- rownames(test_exp_set)
        error_return2 <-
            tryCatch(EWCE::fix_bad_mgi_symbols(
                exp = test_exp_set_char # ,
                # mrk_file_path = sprintf("%s/MRK_List2.rpt", tempdir())
            ),
            error = function(e) e,
            warning = function(w) w
            )

        ## pass data table rather than data frame or matrix
        ## Now detects data.table format and converts to data.frame automatically
        warning_return0 <-
            tryCatch(EWCE::fix_bad_mgi_symbols(
                exp = data.table::data.table(as.data.frame(test_exp_set), keep.rownames = 0)
            ),
            error = function(e) e,
            warning = function(w) w
            )

        # function should warn the user about this -if warning returned function worked
        warning_return <-
            tryCatch(EWCE::fix_bad_mgi_symbols(
                test_exp_set # ,
                # mrk_file_path = sprintf("%s/MRK_List2.rpt", tempdir())
            ),
            error = function(e) e,
            warning = function(w) w
            )

        # running on hgnc rather than mgi should return warnings
        warning_return2 <-
            tryCatch(EWCE::fix_bad_hgnc_symbols(test_exp_set),
                error = function(e) e,
                warning = function(w) w
            )

        # running on hgnc rather than mgi should return warnings
        warning_return3 <-
            tryCatch(
                EWCE::fix_bad_hgnc_symbols(
                    as.data.frame(test_exp_set)
                ),
                error = function(e) e,
                warning = function(w) w
            )

        options(warn = -1)
        hgnc_return <- EWCE::fix_bad_hgnc_symbols(exp = test_exp_set)
        options(warn = 0)

        # Now test if a synonym of a gene in the list is added
        # function should combine them and give sum reads for each sample
        # alt symbol for Tspan12 is Tm4sf12
        rownames(test_exp_set)[7] <- "Tm4sf12"
        EWCE_return <- EWCE::fix_bad_mgi_symbols(
            exp = test_exp_set,
            mrk_file_path = sprintf("%s/MRK_List2.rpt", tempdir()),
            printAllBadSymbols = TRUE
        )
        # Now the returned value should be the average of the two
        sum_exp <- colSums(test_exp_set[
            rownames(test_exp_set) %in% c("Tspan12", "Tm4sf12"),
        ])

        # check nothing changes when there are no issues
        test_exp_set_15 <- test_exp_set[1:5, ]
        EWCE_output_same_input <- EWCE::fix_bad_mgi_symbols(
            test_exp_set_15 # ,
            # mrk_file_path = sprintf("%s/MRK_List2.rpt", tempdir())
        )


        # check input runs on large size of incorrect mgi gene names
        # check it catches warning still but doesn't give an error
        warning_return4 <-
            tryCatch(
                EWCE_output_large <-
                    EWCE::fix_bad_mgi_symbols(
                        cortex_mrna$exp[1:10000, 1:1000] # ,
                        # mrk_file_path = sprintf("%s/MRK_List2.rpt", tempdir())
                    ),
                error = function(e) e,
                warning = function(w) w
            )
        # fail if any of the 3 tests fail
        testthat::expect_true(is(error_return, "error"))
        testthat::expect_true(is(error_return2, "warning"))
        testthat::expect_true(is(warning_return0, "warning"))
        testthat::expect_true(is(warning_return, "warning"))
        testthat::expect_true(is(warning_return2, "warning"))
        testthat::expect_true(is(warning_return3, "warning"))
        testthat::expect_true(EWCE:::is_sparse_matrix(hgnc_return))
        testthat::expect_true(
            all(sum_exp == EWCE_return[rownames(EWCE_return) == "Tspan12", ])
        )
        testthat::expect_true(all.equal(EWCE_output_same_input, test_exp_set_15))
        testthat::expect_true(is(warning_return4, "warning"))

        #----------------------------------------------------------
        # Check SingleCellExperiment gives same results
        # reduce number of genes for speed
        cortex_mrna$exp <- cortex_mrna$exp[seq(1, 6000), ]

        # Make SCE object from SCT
        cortex_mrna_SCE <-
            SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = cortex_mrna$exp),
                colData = cortex_mrna$annot
            )

        # Ensure output is the same for the two data types
        cortex_mrna$exp <- suppressWarnings(
            EWCE::fix_bad_mgi_symbols(exp = cortex_mrna$exp)
        )
        cortex_mrna_SCE <- suppressWarnings(
            EWCE::fix_bad_mgi_symbols(exp = cortex_mrna_SCE)
        )
        # fail if not the exact same
        testthat::expect_true(all.equal(
            SummarizedExperiment::assays(cortex_mrna_SCE)$counts,
            cortex_mrna$exp
        ))
    }
})
