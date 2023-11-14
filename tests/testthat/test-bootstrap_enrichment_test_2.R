test_that("bootstrap enrichment function error handling and geneSizeControl runs as expected", {
    # Test for bootstrap_enrichment_test error handling and geneSizeControl,
    # using sample data in vignette ensure

    if (!is_32bit()) {
        set.seed(12345678)
        # load vignette data
        ctd <- ewceData::ctd()
        human.hits <- example_genelist <- ewceData::example_genelist()
        mouse_to_human_homologs <- ewceData::mouse_to_human_homologs()
        all_hgnc <- ewceData::all_hgnc()

        mouse.hits <- orthogene::convert_orthologs(
            gene_df = human.hits,
            input_species = "human",
            output_species = "mouse",
            gene_output = "dict",
            method = "homologene"
        )
        m2h <- unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])
        mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% example_genelist, "MGI.symbol"])
        mouse.bg <- unique(m2h$MGI.symbol)
        # set input variables
        reps <- 60 # <- Use 60 bootstrap lists so it runs quickly, for publishable analysis use >10000
        level <- 1 # <- Use level 1 annotations (i.e. Interneurons)

        # Test 1: first we want the stop error messages if bad CT controlled passed
        # controlledCT If not NULL, and instead is the name of a cell type,
        # then the bootstrapping controls for expression within that cell type
        controlledCTfail <- "a fake cell type name"

        # should return an error if the function caught the fake cell type
        fail_return <-
            tryCatch(bootstrap_enrichment_test(
                sct_data = ctd,
                hits = unname(mouse.hits),
                genelistSpecies = "mouse",
                sctSpecies = "mouse",
                # bg = mouse.bg,
                reps = reps,
                annotLevel = level,
                controlledCT = controlledCTfail
            ),
            error = function(e) e,
            warning = function(w) w
            )

        # Test 2: second we want to test running bootstrap with geneSizeControl = T
        # geneSizeControl a logical indicating whether you want to control for GC content and
        # transcript length. If set to TRUE then human gene lists should be used rather than mouse.
        # need different dataset from vignette to test
        # human.bg <- unique(c(human.hits, m2h$HGNC.symbol))[1:100] # subset for speed
        results <- EWCE::bootstrap_enrichment_test(
            sct_data = ctd,
            hits = human.hits,
            reps = reps,
            annotLevel = 1,
            geneSizeControl = TRUE,
            genelistSpecies = "human",
            sctSpecies = "mouse",
            output_species = "human"
        )
        # check most significant result
        ewce_sig_cell_types <-
            as.character(results$results[results$results$p <= 0.05 &
                results$results$p == min(results$results$p), "CellType"])

        # Test 3: Next we want to test error handling of check_ewce_genelist_inputs
        # called from bootstrap_enrichment_test
        # test adding a fake mouse hit
        mouse.hits.fake <- c(mouse.hits, "fake")
        fail_return2 <-
            tryCatch(bootstrap_enrichment_test(
                sct_data = ctd,
                hits = mouse.hits.fake,
                genelistSpecies = "mouse",
                sctSpecies = "mouse",
                bg = mouse.bg,
                reps = reps,
                annotLevel = level,
            ),
            error = function(e) e,
            warning = function(w) w
            )

        # test adding a fake mouse bg
        mouse.bg.fake <- mouse.bg[!(mouse.bg %in% all_hgnc)]
        fail_return2_5 <-
            tryCatch(bootstrap_enrichment_test(
                sct_data = ctd,
                hits = mouse.hits,
                bg = mouse.bg.fake,
                reps = reps,
                annotLevel = level,
                genelistSpecies = "mouse",
                sctSpecies = "mouse"
            ),
            error = function(e) e,
            warning = function(w) w
            )

        # Test 4: test adding a sct species other than mouse or human
        rat.hits <- orthogene::convert_orthologs(
            gene_df = human.hits,
            gene_output = "dict",
            input_species = "human",
            output_species = "rat",
            method = "homologene"
        )
        rat_res <- EWCE::bootstrap_enrichment_test(
            sct_data = ctd,
            hits = unname(rat.hits),
            reps = reps,
            annotLevel = level,
            genelistSpecies = "rat",
            sctSpecies = "mouse",
            output_species = "human"
        )
        top_res_rat <- rat_res$results$CellType[1]

        # Test 4: test adding a sct species other than mouse or human
        ctd_monkey <- standardise_ctd(
            ctd = ctd,
            dataset = "ctd_monkey",
            input_species = "mouse",
            output_species = "monkey"
        )
        monkey_res <- bootstrap_enrichment_test(
            sct_data = ctd_monkey,
            hits = mouse.hits,
            reps = reps,
            annotLevel = level,
            sctSpecies = "monkey",
            genelistSpecies = "mouse",
            output_species = "human",
            store_gene_data = FALSE
        )
        top_res_monkey <- monkey_res$results$CellType[1]


        # Test 4.5: test adding a fake species that doesn't exist in orthogene
        fail_return3 <-
            tryCatch(bootstrap_enrichment_test(
                sct_data = ctd,
                hits = mouse.hits,
                reps = reps,
                annotLevel = level,
                genelistSpecies = "mouse",
                sctSpecies = "godzilla",
                output_species = "human"
            ),
            error = function(e) e,
            warning = function(w) w
            )

        # Test 5: test supplying less than 4 hits for mouse
        # fail_return4 <-
        #    tryCatch(bootstrap_enrichment_test(
        #        sct_data = ctd, hits = mouse.hits[1:3],
        #        bg = mouse.bg, reps = reps, annotLevel = level
        #    ),
        #    error = function(e) e,
        #    warning = function(w) w
        #    )
        # fail_return4 REMOVED TO SAVE ON RUN TIME OF CHECK

        # Test 6: test supplying less than 4 hits for human
        fail_return5 <-
            tryCatch(bootstrap_enrichment_test(
                sct_data = ctd,
                hits = mouse.hits[seq(1, 3)],
                reps = reps,
                annotLevel = level,
                genelistSpecies = "human",
                sctSpecies = "mouse",
                output_species = "human"
            ),
            error = function(e) e,
            warning = function(w) w
            )

        # test 7: running human species with mouse gene list
        # should fail because of gene names
        fail_return6 <-
            tryCatch(bootstrap_enrichment_test(
                sct_data = ctd,
                hits = mouse.hits,
                bg = mouse.bg,
                reps = reps,
                annotLevel = 1,
                geneSizeControl = TRUE,
                genelistSpecies = "mouse",
                sctSpecies = "human"
            ),
            error = function(e) e,
            warning = function(w) w
            )
        # ctd_fake <- ctd
        # just rename genes to hgnc names
        # rownames(ctd_fake[[1]]$mean_exp) <-
        #    all_hgnc[seq_len(length(rownames(ctd[[1]]$mean_exp)))]
        # should fail as geneSizeControl set to true
        # fail_return7 <-
        #    tryCatch(bootstrap_enrichment_test(
        #        sct_data = ctd_fake, hits = mouse.hits,
        #        bg = mouse.bg, reps = reps, annotLevel = 1,
        #        geneSizeControl = TRUE, genelistSpecies = "mouse",
        #        sctSpecies = "human"
        #    ),
        #    error = function(e) e,
        #    warning = function(w) w
        #    )
        # fail_return7 REMOVED TO SAVE ON RUN TIME OF CHECK
        # fail if any subtest fails

        #### Run tests ####
        testthat::expect_true(methods::is(fail_return, "error"))
        testthat::expect_true("microglia" %in% ewce_sig_cell_types)
        testthat::expect_true(methods::is(fail_return2, "error"))
        testthat::expect_true(top_res_rat == "microglia")
        testthat::expect_false(any(is.na(rat_res$results)))
        testthat::expect_true(top_res_monkey == "microglia")
        testthat::expect_false(any(is.na(monkey_res$results)))
        testthat::expect_true(methods::is(fail_return3, "error"))
        testthat::expect_true(methods::is(fail_return2_5, "error"))
        testthat::expect_true(methods::is(fail_return5, "error"))
        testthat::expect_true(methods::is(fail_return6, "error"))
    }
})
