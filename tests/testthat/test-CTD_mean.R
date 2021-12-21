test_that("Correct mean expression values of CTD calculated", {

    # Test for mean calculation from calculate_meanexp_for_level
    # in generate_celltype_data

    if (!is_32bit()) {
        # create some fake exp data
        set.seed(101)
        rand_data <- runif(1000)
        exp_set <- matrix(ifelse(rand_data < 0.55, 0, rand_data),
            ncol = 10,
            dimnames = list(
                paste0("gene_", seq(1, 100)),
                paste0("sample_", seq(1, 10))
            )
        )
        # Make exp into a sparse matrix
        exp_set <- Matrix::Matrix(exp_set)
        # now create 4 groups for the ten samples to get the means of
        grouping <- vector(mode = "list", length = 4)
        # just creating four levels so the test can be run four times
        names(grouping) <- paste0("level", seq(1, 4), "class")
        grouping <- lapply(grouping, function(x) {
            as.character(sample(paste0("celltype_", seq(1, 4)),
                10,
                replace = TRUE
            ))
        })

        # Now apply EWCE function
        fNames <- EWCE::generate_celltype_data(
            exp = exp_set,
            annotLevels = grouping,
            ## IMPORTANT!: Must be false so fakes genes don't all get dropped.
            convert_orths = FALSE,
            groupName = "testthat",
            savePath = tempdir(),
            ### IMPORTANT! Testthat causes conflicts with parallelization
            no_cores = 1
        )
        # load res - named ctd
        ctd <- EWCE::load_rdata(fNames)
        EWCE_res_meanexp <- lapply(ctd, function(x) {
            methods::as(x$mean_exp, "sparseMatrix")
        })

        # check alternative method
        alt_res <- vector(
            mode = "list",
            length = length(ctd)
        )
        for (i in seq_len(length(grouping))) {
            group_i <- grouping[[i]]
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
            alt_res[[i]] <- methods::as(unique_cell_type_means, "sparseMatrix")
            # fail if mean expression values aren't the same across the 3 tests
            testthat::expect_true(all.equal(
                EWCE_res_meanexp[[i]],
                alt_res[[i]]
            ))
        }
    }
})
