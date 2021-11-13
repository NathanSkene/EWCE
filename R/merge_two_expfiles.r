#' Merge two exp files
#'
#' \code{merge_two_expfiles} Used to combine two single cell type datasets.
#'
#' @param exp1 Numerical expression matrix for dataset1 with row for each gene
#' and column for each cell. Row names are gene symbols. Column names
#' are cell IDs which can be cross referenced against the annot data frame.
#' @param exp2 Numerical expression matrix for dataset2 with row for each gene
#' and column for each cell. Row names are gene symbols. Column names
#' are cell IDs which can be cross referenced against the annot data frame.
#' @param annot1 Annotation data frame for dataset1 which contains three
#' columns at least: cell_id, level1class and level2class
#' @param annot2 Annotation data frame for dataset2 which contains three
#' columns at least: cell_id, level1class and level2class
#' @param name1 Name used to refer to dataset 1. Leave blank if it's already a
#' merged dataset.
#' @param name2 Name used to refer to dataset 2. Leave blank if it's already a
#' merged dataset.
#' @param as_sparse Convert the merged \code{exp} to a sparse matrix.
#' @param as_DelayedArray Convert the merged \code{exp} to
#' a \code{DelayedArray}.
#' @param verbose Print messages.
#'
#' @return List containing merged exp and annot.
#'
#' @examples
#' cortex_mrna <- ewceData::cortex_mrna()
#' exp1 <- cortex_mrna$exp[, 1:50]
#' exp2 <- cortex_mrna$exp[, 51:100]
#' annot1 <- cortex_mrna$annot[1:50, ]
#' annot2 <- cortex_mrna$annot[51:100, ]
#' merged_res <- EWCE::merge_two_expfiles(
#'     exp1 = exp1,
#'     exp2 = exp2,
#'     annot1 = annot1,
#'     annot2 = annot2,
#'     name1 = "dataset1",
#'     name2 = "dataset2"
#' )
#' @export
merge_two_expfiles <- function(exp1,
                               exp2,
                               annot1,
                               annot2,
                               name1 = "",
                               name2 = "",
                               as_sparse = TRUE,
                               as_DelayedArray = FALSE,
                               verbose = TRUE) {
    all_genes <- intersect(unique(rownames(exp1)), unique(rownames(exp2)))
    removed_genes1 <- rownames(exp1)[!rownames(exp1) %in% rownames(exp2)]
    removed_genes2 <- rownames(exp2)[!rownames(exp2) %in% rownames(exp1)]
    messager(formatC(length(removed_genes1), big.mark = ","),
        "non-overlapping gene(s) removed from exp1.",
        v = verbose
    )
    messager(formatC(length(removed_genes2), big.mark = ","),
        "non-overlapping gene(s) removed from exp2.",
        v = verbose
    )
    messager(formatC(length(all_genes), big.mark = ","),
        "intersecting genes remain.",
        v = verbose
    )

    #### Check for correct column headers in annot data.frames ####
    if (sum(c("cellid") %in% colnames(annot1)) == 1) {
        colnames(annot1)[colnames(annot1) == "cellid"] <- "cell_id"
    }
    if (sum(c("cellid") %in% colnames(annot2)) == 1) {
        colnames(annot2)[colnames(annot2) == "cellid"] <- "cell_id"
    }
    err_msg <- paste0(
        "ERROR: annot1 doesn't have either cell_id,",
        " level1class or level2class columns"
    )
    err_msg2 <- paste0(
        "ERROR: annot1 doesn't have either cell_id,",
        " level1class or level2class columns"
    )
    if (!sum(c("cell_id", "level1class", "level2class") %in%
        colnames(annot1)) == 3) {
        stop(err_msg)
    }
    if (!sum(c("cell_id", "level1class", "level2class") %in%
        colnames(annot2)) == 3) {
        stop(err_msg2)
    }

    # If one of the exp matrices is really a matrix (not a dataframe)
    # the code won't work... so force conversion
    exp1 <- to_dataframe(X = exp1)
    exp2 <- to_dataframe(X = exp2)
    #### Merge the expression dfs, setting undetected genes to 0 ####
    exp1b <- exp1[all_genes, ]
    exp1b[is.na(exp1b)] <- 0
    rownames(exp1b) <- all_genes
    exp2b <- exp2[all_genes, ]
    exp2b[is.na(exp2b)] <- 0
    rownames(exp2b) <- all_genes
    exp <- cbind(exp1b, exp2b)
    remove(exp1, exp2, exp1b, exp2b)

    # Ensure important annotation columns are not stored as factors
    annot1$cell_id <- as.character(annot1$cell_id)
    annot2$cell_id <- as.character(annot2$cell_id)
    annot1$level1class <- as.character(annot1$level1class)
    annot1$level2class <- as.character(annot1$level2class)
    if (is.null(annot1$dataset_name)) {
        annot1$dataset_name <- name1
    }
    annot2$level1class <- as.character(annot2$level1class)
    annot2$level2class <- as.character(annot2$level2class)
    if (is.null(annot2$dataset_name)) {
        annot2$dataset_name <- name2
    }
    keepTISSUE <- FALSE
    if (("tissue" %in% colnames(annot1)) & ("tissue" %in% colnames(annot2))) {
        keepTISSUE <- TRUE
        annot1$tissue <- as.character(annot1$tissue)
        annot2$tissue <- as.character(annot2$tissue)
    }

    # Setup new annotation dataframe
    cell_id <- c(annot1$cell_id, annot2$cell_id)
    level1class <- c(annot1$level1class, annot2$level1class)
    level2class <- c(annot1$level2class, annot2$level2class)
    dataset_name <- c(
        as.character(annot1$dataset_name),
        as.character(annot2$dataset_name)
    )
    if (!keepTISSUE) {
        annot <- data.frame(
            cell_id = cell_id, level1class = level1class,
            level2class = level2class,
            dataset_name = dataset_name
        )
    } else {
        tissue <- c(annot1$tissue, annot2$tissue)
        annot <- data.frame(
            cell_id = cell_id, level1class = level1class,
            level2class = level2class,
            dataset_name = dataset_name, tissue = tissue
        )
    }


    # Drop expression data without annotation data
    numMissingAnnot <- dim(exp)[2] - sum(annot$cell_id %in% colnames(exp))
    if (numMissingAnnot > 0) {
        txt <- paste0(
            "Warning: %s cells are missing annotation data",
            " and have been dropped"
        )
        message(sprintf(txt, numMissingAnnot))
        exp <- exp[, as.character(annot$cell_id)]
    }
    #### Convert to sparse matrix (optional) ####
    exp <- to_sparse_matrix(
        exp = exp,
        as_sparse = as_sparse,
        verbose = verbose
    )
    #### Convert to DelayedArray (optional) ####
    exp <- to_delayed_array(
        exp = exp,
        as_DelayedArray = as_DelayedArray,
        verbose = verbose
    )
    return(list(
        exp = exp,
        annot = annot
    ))
}
