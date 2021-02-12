# Test for bootstrap.enrichment.test error handling and geneSizeControl, using sample data in vignette ensure
test_that("bootstrap enrichment function error handling and geneSizeControl runs as expected", {
    # load vignette data
    data("ctd")
    data(example_genelist)
    data("mouse_to_human_homologs")
    m2h <- unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])
    mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% example_genelist, "MGI.symbol"])
    # mouse.bg  = unique(setdiff(m2h$MGI.symbol,mouse.hits))
    mouse.bg <- unique(m2h$MGI.symbol)
    # set input variables
    reps <- 100 # <- Use 100 bootstrap lists so it runs quickly, for publishable analysis use >10000
    level <- 1 # <- Use level 1 annotations (i.e. Interneurons)

    # Test 1: first we want the stop error messages if bad CT controlled passed
    # controlledCT If not NULL, and instead is the name of a cell type,
    # then the bootstrapping controls for expression within that celltype
    controlledCTfail <- "a fake celltype name"

    # should return an error if the function caught the fake cell type
    fail_return <-
        tryCatch(bootstrap.enrichment.test(
            sct_data = ctd, hits = mouse.hits,
            bg = mouse.bg, reps = reps, annotLevel = level,
            controlledCT = controlledCTfail
        ),
        error = function(e) e,
        warning = function(w) w
        )

    # Test 2: second we want to test running bootstrap with geneSizeControl = T
    # geneSizeControl a logical indicating whether you want to control for GC content and
    # transcript length. If set to TRUE then human gene lists should be used rather than mouse.
    # need different dataset from vignette to test
    human.hits <- example_genelist
    human.bg <- unique(c(human.hits, m2h$HGNC.symbol))
    set.seed(12345678)
    results <- bootstrap.enrichment.test(
        sct_data = ctd, hits = human.hits,
        bg = human.bg, reps = reps, annotLevel = 1,
        geneSizeControl = TRUE, genelistSpecies = "human",
        sctSpecies = "mouse"
    )
    # check most significant result
    ewce_sig_cell_types <-
        as.character(results$results[results$results$p <= 0.05 &
            results$results$p == min(results$results$p), "CellType"])

    # Test 3: Next we want to test error handling of check.ewce.genelist.inputs called from bootstrap.enrichment.test
    # test adding a fake mouse hit
    mouse.hits.fake <- c(mouse.hits, "fake")
    fail_return2 <-
        tryCatch(bootstrap.enrichment.test(
            sct_data = ctd, hits = mouse.hits.fake,
            bg = mouse.bg, reps = reps, annotLevel = level,
        ),
        error = function(e) e,
        warning = function(w) w
        )

    # test adding a fake mouse bg
    mouse.bg.fake <- mouse.bg[!(mouse.bg %in% EWCE::all_hgnc)]
    fail_return2_5 <-
        tryCatch(bootstrap.enrichment.test(
            sct_data = ctd, hits = mouse.hits,
            bg = mouse.bg.fake, reps = reps, annotLevel = level,
            genelistSpecies = "human"
        ),
        error = function(e) e,
        warning = function(w) w
        )

    # Test 4: test adding a sct species other than mouse or human
    fail_return3 <-
        tryCatch(bootstrap.enrichment.test(
            sct_data = ctd, hits = mouse.hits,
            bg = mouse.bg, reps = reps, annotLevel = level,
            sctSpecies = "rat"
        ),
        error = function(e) e,
        warning = function(w) w
        )

    # Test 5: test supplying less than 4 hits for mouse
    fail_return4 <-
        tryCatch(bootstrap.enrichment.test(
            sct_data = ctd, hits = mouse.hits[1:3],
            bg = mouse.bg, reps = reps, annotLevel = level
        ),
        error = function(e) e,
        warning = function(w) w
        )

    # Test 6: test supplying less than 4 hits for human
    fail_return5 <-
        tryCatch(bootstrap.enrichment.test(
            sct_data = ctd, hits = mouse.hits[1:3],
            bg = mouse.bg, reps = reps, annotLevel = level,
            genelistSpecies = "human"
        ),
        error = function(e) e,
        warning = function(w) w
        )

    # test 7: running human species with mouse gene list
    # should fail because of gene names
    fail_return6 <-
        tryCatch(bootstrap.enrichment.test(
            sct_data = ctd, hits = mouse.hits,
            bg = mouse.bg, reps = reps, annotLevel = 1,
            geneSizeControl = TRUE, genelistSpecies = "mouse",
            sctSpecies = "human"
        ),
        error = function(e) e,
        warning = function(w) w
        )
    ctd_fake <- ctd
    # just rename genes to hgnc names
    rownames(ctd_fake[[1]]$mean_exp) <- EWCE::all_hgnc[seq_len(length(rownames(ctd[[1]]$mean_exp)))]
    # should fail as geneSizeControl set to true
    fail_return7 <-
        tryCatch(bootstrap.enrichment.test(
            sct_data = ctd_fake, hits = mouse.hits,
            bg = mouse.bg, reps = reps, annotLevel = 1,
            geneSizeControl = TRUE, genelistSpecies = "mouse",
            sctSpecies = "human"
        ),
        error = function(e) e,
        warning = function(w) w
        )
    # fail if any subtest fails
    expect_equal(all(
        # Test 1
        is(fail_return, "error"),
        # Test 2
        ewce_sig_cell_types == "microglia",
        # Test 3
        is(fail_return2, "error"),
        # Test 3
        is(fail_return2_5, "error"),
        # Test 4
        is(fail_return3, "error"),
        # Test 5
        is(fail_return4, "error"),
        # Test 6
        is(fail_return5, "error"),
        # test 7
        is(fail_return6, "error"),
        is(fail_return7, "error")
    ), TRUE)
})
