# Test for mean calculation from calculate.meanexp.for.level in generate.celltype.data
test_that("Correct mean expression values of CDT calculated", {
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
        groupName = "testthat"
    )
    # load res - named ctd
    load(fNames)
    EWCE_res_meanexp <- lapply(ctd, function(x) x$mean_exp)

    # check alternative method
    alt_res <- vector(mode = "list", length = 3)
    i <- 1
    for (group_i in grouping) {
        unique_cell_types <- unique(group_i)
        unique_cell_type_means <- matrix(
            ncol = 0,
            nrow = nrow(exp_set),
            dimnames = list(rownames(exp_set))
        )
        for (uct_i in sort(unique_cell_types)) {
            unique_cell_type_means <- cbind(
                unique_cell_type_means,
                rowSums(as.matrix(exp_set[, group_i == uct_i])) /
                    sum(group_i == uct_i)
            )
        }
        colnames(unique_cell_type_means) <- sort(unique_cell_types)
        alt_res[[i]] <- unique_cell_type_means
        i <- i + 1
    }

    # fail if mean expression values aren't the same across the 3 tests
    expect_equal(all.equal(EWCE_res_meanexp, alt_res), TRUE)
})
