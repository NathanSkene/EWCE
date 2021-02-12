# test ewce_expression_data which calls bootstrap.enrichment.test
# Test for specificity return from generate.celltype.data
test_that("EWCE expression data creation", {
    # bootstrap.enrichment.test tested in depth separately, testing returned specificity here
    set.seed(12345678)
    data(tt_alzh)
    tt_results <- ewce_expression_data(sct_data = ctd, tt = tt_alzh, annotLevel = 1, ttSpecies = "human", sctSpecies = "mouse")
    up_celltypes <- c(
        names(tt_results$hit.cells.up)[tt_results$hit.cells.up == max(tt_results$hit.cells.up)],
        names(tt_results$hit.cells.up)[tt_results$hit.cells.up == min(tt_results$hit.cells.up)]
    )
    dwn_celltypes <- c(
        names(tt_results$hit.cells.down)[tt_results$hit.cells.down == max(tt_results$hit.cells.down)],
        names(tt_results$hit.cells.down)[tt_results$hit.cells.down == min(tt_results$hit.cells.down)]
    )

    # check that max and min specificity matches for up and down
    known_maxmin_up_celltypes <- c("endothelial-mural", "interneurons")
    known_maxmin_dwn_celltypes <- c("pyramidal SS", "oligodendrocytes")

    # fail if specificity max and min cell types doesn't match
    expect_equal(
        (all(up_celltypes == known_maxmin_up_celltypes) &&
            all(dwn_celltypes == known_maxmin_dwn_celltypes)),
        TRUE
    )
})
