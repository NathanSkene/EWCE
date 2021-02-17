#' generate.celltype.data
#'
#' \code{generate.celltype.data} Takes expression & cell type annotations and 
#' creates celltype_data files which contain the mean and specificity matrices
#'
#' @param exp Numerical matrix with row for each gene and column for each cell. 
#' Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be 
#' cross referenced against the annot data frame.
#' @param annotLevels List with arrays of strings containing the cell type 
#' names associated with each column in exp
#' @param groupName A human readable name for refering to the dataset being 
#' loaded
#' @param no_cores Number of cores that should be used to speedup the 
#' computation. Use no_cores = 1 when using this package in windows system.
#' @param savePath Directory where the CTD file should be saved
#' @param normSpec Boolean indicating whether specificity data should be 
#' transformed to a normal distribution by cell type, giving equivalent scores 
#' across all cell types.
#' @return Filenames for the saved celltype_data files
#' @examples
#' # Load the single cell data
#' data("cortex_mrna", package="ewceData")
#' expData <- cortex_mrna$exp
#' expData <- expData[1:500, ] # Use only a subset to keep the example quick
#' l1 <- cortex_mrna$annot$level1class
#' l2 <- cortex_mrna$annot$level2class
#' annotLevels <- list(l1 = l1, l2 = l2)
#' fNames_ALLCELLS <- 
#'     generate.celltype.data(exp = expData, annotLevels, "allKImouse")
#' @export
#' @import parallel
#' @import future
#' @import ggdendro
#' @import gridExtra
#' @importFrom Matrix Matrix
#' @import RNOmni
#' @import ggdendro

generate.celltype.data <- function(exp, annotLevels, groupName, no_cores = 1, 
                                    savePath = "~/", normSpec = FALSE) {
    err_msg <- paste0("NA values detected in expresson matrix. All NA values",
                        " should be removed before calling EWCE.")
    if (sum(is.na(exp)) > 0) {
        stop(err_msg)
    }
    # Calculate the number of cores
    # cl <- parallel::makeCluster(no_cores)
    # print(sprintf("Using %s cores",no_cores))
    # First, check the number of annotations equals the number of columns 
    # in the expression data
    err_msg2 <- paste0("Error: length of all annotation levels must equal",
                        " the number of columns in exp matrix")
    lapply(annotLevels, test <- function(x, exp) {
        if (length(x) != dim(exp)[2]) {
            stop(err_msg2)
        }
    }, exp)

    # Check group name
    err_msg3 <- paste0("ERROR: groupName must be set. groupName is used to",
                        " label the files created by this function.")
    if (is.null(groupName)) {
        stop(err_msg3)
    }
    if (groupName == "") {
        stop(err_msg3)
    }

    # Convert characters to numbers
    if (!is(exp)[1] == "dgCMatrix") {
        exp <- suppressWarnings(apply(exp, 2, function(x) {
            storage.mode(x) <- "double"
            x
        }))
    }

    # Make exp into a sparse matrix
    exp <- Matrix::Matrix(exp)
    ctd <- list()
    for (i in seq_len(length(annotLevels))) {
        ctd[[length(ctd) + 1]] <- list(annot = annotLevels[[i]])
    }
    ctd2 <- mclapply(ctd, calculate.meanexp.for.level, exp, mc.cores = no_cores)
    ctd3 <- mclapply(ctd2, calculate.specificity.for.level, mc.cores = no_cores)
    ctd <- ctd3
    # stopCluster(cl)

    # apply rank norm transformation if hyp param set
    if (isTRUE(normSpec)) {ctd <- lapply(ctd, rNorm)}
    # ADD DENDROGRAM DATA TO CTD
    ctd <- lapply(ctd, bin.specificity.into.quantiles, numberOfBins = 40)
    ctd <- lapply(ctd, prep.dendro)

    fNames <- sprintf("%s/CellTypeData_%s.rda", savePath, groupName)
    save(ctd, file = fNames)
    return(fNames)
}

#Helper functions
calculate.meanexp.for.level <- function(ctd_oneLevel, expMatrix) {
    err_msg <- paste0("There are an equal number of celltypes in expMatrix",
                        " and ctd_oneLevel but the names do not match")
    if (dim(expMatrix)[2] == length(unique(ctd_oneLevel$annot))) {
        print(dim(expMatrix)[2])
        print(length(ctd_oneLevel$annot))
        if (sum(!colnames(expMatrix) == ctd_oneLevel$annot) != 0) {
            stop(err_msg)
        }
        ctd_oneLevel$mean_exp <- as.matrix(expMatrix)
    } else {
        # Sum reads in each cell type
        # mean_exp = apply(expMatrix,1,aggregate.over.celltypes,
        #                   ctd_oneLevel$annot)
        mm <- model.matrix(~ 0 + ctd_oneLevel$annot)
        colnames(mm) <- names(table(ctd_oneLevel$annot))
        mat.summary.mm1 <- expMatrix %*% mm
        
        # Divide by the number of cells to get the mean
        cellCounts <- table(ctd_oneLevel$annot)
        for (i in seq_len(dim(mat.summary.mm1)[2])) {
            mat.summary.mm1[, i] <- mat.summary.mm1[, i] / cellCounts[i]
        }
        
        ctd_oneLevel$mean_exp <- as.matrix(mat.summary.mm1)
    }
    return(ctd_oneLevel)
}

calculate.specificity.for.level <- function(ctd_oneLevel) {
    normalised_meanExp <- t(t(ctd_oneLevel$mean_exp) * 
                                (1 / colSums(ctd_oneLevel$mean_exp)))
    ctd_oneLevel$specificity <- normalised_meanExp / 
        (apply(normalised_meanExp, 1, sum) + 0.000000000001)
    return(ctd_oneLevel)
}

# Use the rank norm transformation on specificity
rNorm <- function(ctdIN) {
    ctdIN$specificity <- apply(ctdIN$specificity, 2, RNOmni::RankNorm)
    return(ctdIN)
}