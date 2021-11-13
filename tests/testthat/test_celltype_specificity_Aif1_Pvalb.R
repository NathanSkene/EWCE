test_that("Cell type specificity Aif1, Pvalb", {
    if (!is_32bit()) {
        # Test Cell type specificity calculations Aif1, Pvalb
        # inspect specificity matrix for type 1 cell types from
        # generate_celltype_data
        # Use Vignette Dataset to check function output
        cortex_mrna <- ewceData::cortex_mrna()
        nKeep <- 1000
        must_keep <- c("Apoe", "Gfap", "Gapdh", "Pvalb", "Aif1")
        set.seed(123458)
        keep_genes <- c(must_keep, sample(rownames(cortex_mrna$exp), 995))
        cortex_mrna$exp <- cortex_mrna$exp[keep_genes, ]

        # Note: not normalising data to save time
        exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(
            exp = cortex_mrna$exp, # _scT_normed,
            level2annot = cortex_mrna$annot$level2class,
            drop_nonhuman_genes = FALSE,
            ### IMPORTANT! Testthat causes conflicts with parallelization
            no_cores = 1
        )
        annotLevels <- list(
            level1class = cortex_mrna$annot$level1class,
            level2class = cortex_mrna$annot$level2class
        )
        fNames_CortexOnly <- EWCE::generate_celltype_data(
            exp = exp_CortexOnly_DROPPED,
            annotLevels = annotLevels,
            convert_orths = FALSE, ### Keep original mouse genes
            input_species = "mouse",
            output_species = "human",
            groupName = "kiCortexOnly",
            savePath = tempdir(),
            ### IMPORTANT! Testthat causes conflicts with parallelization
            no_cores = 1
        )
        # load and inspect specificity matrix
        ctd <- EWCE::load_rdata(fNames_CortexOnly[1])
        #### Use mouse orthologs to query
        EWCE_return_Pvalb <- sort(ctd[[1]]$specificity["Pvalb", ],
            decreasing = TRUE
        )
        EWCE_return_Aif1 <- sort(ctd[[1]]$specificity["Aif1", ],
            decreasing = TRUE
        )

        pvalb_interneurons <- names(EWCE_return_Pvalb)[[1]] == "interneurons"
        aif1_microglia <- names(EWCE_return_Aif1)[[1]] == "microglia"

        # Try running with rank norm transformation normSpec set to True
        # Ensure same most specific cell types remain
        ctd_not_norm_spec <- ctd
        fNames_CortexOnly <- EWCE::generate_celltype_data(
            exp = exp_CortexOnly_DROPPED,
            annotLevels = annotLevels,
            convert_orths = FALSE, ### Keep original mouse genes
            input_species = "mouse",
            output_species = "human",
            groupName = "kiCortexOnly",
            normSpec = TRUE,
            savePath = tempdir(),
            ### IMPORTANT! Testthat causes conflicts with parallelization
            no_cores = 1
        )
        # load and inspect specificity matrix
        ctd <- EWCE::load_rdata(fNames_CortexOnly[1])

        EWCE_return_Pvalb_norm <- sort(ctd[[1]]$specificity["Pvalb", ],
            decreasing = TRUE
        )
        EWCE_return_Aif1_norm <- sort(ctd[[1]]$specificity["Aif1", ],
            decreasing = TRUE
        )

        pvalb_interneurons_norm <-
            names(EWCE_return_Pvalb_norm)[[1]] == "interneurons"
        aif1_microglia_norm <-
            names(EWCE_return_Aif1_norm)[[1]] == "microglia"

        # also check the results appear the same
        # - apart from the transformed, normalised specificity values
        ctd_comparable <- all(
            all.equal(ctd_not_norm_spec[[1]]$annot, ctd[[1]]$annot),
            all.equal(ctd_not_norm_spec[[1]]$mean_exp, ctd[[1]]$mean_exp),
            ### orthogene stores convert genes as named list.
            # Must unname to make comparable.
            all.equal(
                unname(rownames(ctd_not_norm_spec[[1]]$specificity)),
                rownames(ctd[[1]]$specificity)
            ),
            all.equal(
                colnames(ctd_not_norm_spec[[1]]$specificity),
                colnames(ctd[[1]]$specificity)
            )
        )
        # Fail if most specific cell type for
        # Aif1, Pvalb not microglia and interneurons respectively
        testthat::expect_true(pvalb_interneurons)
        testthat::expect_true(aif1_microglia)
        testthat::expect_true(pvalb_interneurons_norm)
        testthat::expect_true(aif1_microglia_norm)
        testthat::expect_true(ctd_comparable)
    }
})
