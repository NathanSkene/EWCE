# Test for bootstrap.enrichment.test, using sample data in vignette ensure
# gives microglia as the only significant enrichment
test_that("bootstrap enrichment runs as expected", {
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
    set.seed(12345678)
    full_results <-
        bootstrap.enrichment.test(
            sct_data = ctd, hits = mouse.hits,
            bg = mouse.bg, reps = reps, annotLevel = level
        )

    ewce_sig_cell_types <- as.character(full_results$results[full_results$results$p < 0.05, "CellType"])
    # fail if any but microglia returned i.e. nothing or multiples too
    expect_equal(ewce_sig_cell_types, "microglia")
})
