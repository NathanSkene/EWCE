# Test for specificity return from generate.celltype.data
test_that("Correct specificity values of CDT calculated", {
    # create some fake exp data
    set.seed(101)
    rand_data <- runif(1000)
    exp_set <- matrix(ifelse(rand_data < 0.55, 0, rand_data),
        ncol = 10,
        dimnames = list(paste0("gene_", 1:100), paste0("sample_", 1:10))
    )
    # Make exp into a sparse matrix
    exp_set <- Matrix::Matrix(exp_set)
    # now create 4 groups for the ten samples to get the means of
    grouping <- vector(mode = "list", length = 4)
    # just creating four levels so the test can be run four times
    names(grouping) <- paste0("level", 1:4, "class")
    grouping <- lapply(grouping, function(x) {
          as.character(sample(paste0("celltype_", 1:4), 10, replace = TRUE))
      })

    # Now apply EWCE function
    fNames <- EWCE::generate.celltype.data(
        exp = exp_set,
        annotLevels = grouping,
        groupName = "testthat",
        savePath = tempdir()
    )
    # load res - named ctd
    load(fNames)
    EWCE_res_spec <- lapply(ctd, function(x) x$specificity)

    # check against known answer for gene 16
    test <- setNames(c(1, 2), c("A", "B"))
    known_ans <- list(
        setNames(
            c(0.4450643, 0.2458045, 0.3091312),
            c("celltype_2", "celltype_3", "celltype_4")
        ),
        setNames(
            c(0.4050525, 0.4271945, 0.1677530),
            c("celltype_1", "celltype_2", "celltype_4")
        ),
        setNames(
            c(0.2022591, 0.3428756, 0.1924686, 0.2623967),
            c("celltype_1", "celltype_2", "celltype_3", "celltype_4")
        ),
        setNames(
            c(0.1583313, 0.2530405, 0.3367327, 0.2518955),
            c("celltype_1", "celltype_2", "celltype_3", "celltype_4")
        )
    )
    EWCE_res_gene_16 <- lapply(EWCE_res_spec, function(x) round(x[rownames(x) == "gene_16", ], 7))

    # fail if specificity values aren't the same for gene 16 across the 3 tests
    expect_equal(all.equal(EWCE_res_gene_16, known_ans), TRUE)
})
