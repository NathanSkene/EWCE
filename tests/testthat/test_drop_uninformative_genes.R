# Test for add.res.to.merging.list, using sample data in vignette
test_that("Drop uninformative genes", {
    # inspect drop.uninformative.genes
    # Use Vignette Dataset to check function output
    cortex_mrna <- cortex_mrna()
    if (!file.exists("MRK_List2.rpt")) {
        download.file("http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt", destfile = "MRK_List2.rpt")
    }

    nKeep <- 1000
    must_keep <- c("Apoe", "Gfap", "Gapdh")
    set.seed(123458)
    keep_genes <- c(must_keep, sample(rownames(cortex_mrna$exp), 997))
    cortex_mrna$exp <- cortex_mrna$exp[keep_genes, ]

    # NOTE not normalising for this test, too time intensive

    exp_CortexOnly_DROPPED <- drop.uninformative.genes(
        exp = cortex_mrna$exp,
        level2annot = cortex_mrna$annot$level2class
    )

    # remove folder once tested
    unlink("MRK_List2.rpt", recursive = TRUE)

    # check the number of dropped genes is the same as expected from past runs: 248 else fail
    expect_equal(nrow(cortex_mrna$exp) - nrow(exp_CortexOnly_DROPPED), 248)
})
