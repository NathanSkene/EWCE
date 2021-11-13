#' fix_bad_hgnc_symbols
#'
#' Given an expression matrix, wherein the rows are supposed to be HGNC
#' symbols, find those symbols which are not official HGNC symbols, then
#' correct them if possible. Return the expression matrix with corrected
#' symbols.
#'
#' @param exp An expression matrix where the rows are HGNC symbols or a
#' SingleCellExperiment (SCE) or other
#'  Ranged Summarized Experiment (SE) type object.
#' @param dropNonHGNC Boolean. Should symbols not recognised as HGNC symbols
#' be dropped?
#' @param as_sparse Convert \code{exp} to sparse matrix.
#' @param verbose Print messages.
#'
#' @returns Returns the expression matrix with the rownames corrected and rows
#' representing the same gene merged. If a SingleCellExperiment (SCE) or other
#' Ranged Summarized Experiment (SE) type object was inputted this will be
#' returned with the corrected expression matrix under counts.
#'
#' @examples
#' # create example expression matrix, could be part of a exp, annot list obj
#' exp <- matrix(data = runif(70), ncol = 10)
#' # Add HGNC gene names but add with an error:
#' # MARCH8 is a HGNC symbol which if opened in excel will convert to Mar-08
#' rownames(exp) <-
#'     c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1", "MT-ND1", "Mar-08")
#' exp <- fix_bad_hgnc_symbols(exp)
#' # fix_bad_hgnc_symbols warns the user of this possible issue
#' @export
#' @import HGNChelper
#' @import ewceData
#' @importFrom SummarizedExperiment rowRanges assays
#' @importFrom methods as
fix_bad_hgnc_symbols <- function(exp,
                                 dropNonHGNC = FALSE,
                                 as_sparse = TRUE,
                                 verbose = TRUE) {
    # Check arguments
    err_msg <- paste0(
        "'exp' is null. It should be a numerical",
        " matrix with the rownames being HGNC symbols."
    )
    if (is.null(exp)) {
        stop(err_msg)
    }
    # Check if input is an SCE or SE object
    res_sce <- check_sce(exp)
    exp <- res_sce$exp
    SE_obj <- res_sce$SE_obj
    #### Convert to data.table --> data.frame ####
    exp <- dt_to_df(exp = exp)
    #### Convert characters to numbers ####
    exp <- check_numeric(exp = exp)
    # Check that exp is not some weird input format like those generated
    # by readr functions
    if (!any(
        is_sparse_matrix(exp),
        is_matrix(exp),
        is_delayed_array(exp),
        is.data.frame(exp)
    )) {
        stop("exp must be either matrix (sparse or dense) or data.frame")
    }
    #### First, find which gene symbols are not proper HGNC symbols ####
    all_hgnc <- ewceData::all_hgnc()
    not_HGNC <- rownames(exp)[!rownames(exp) %in% all_hgnc]
    messager(sprintf(
        "%s of %s are not proper HGNC symbols.",
        length(unique(not_HGNC)), dim(exp)[1]
    ), v = verbose)

    #### If dropNonHGNC==TRUE then dropNonHGNC symbols ####
    if (dropNonHGNC == TRUE) {
        exp <- exp[!rownames(exp) %in% not_HGNC, ]
    }
    #### First, check if any gene symbols have been corrupted by excel ####
    date_like <- not_HGNC[grep("SEP|MAR|FEB|DEC|Sep|Mar|Feb|Dec", not_HGNC)]
    if (length(date_like) > 0) {
        messager(sprintf(
            "Possible corruption of gene names by excel: %s",
            paste(date_like, collapse = ", ")
        ), v = verbose)
        warning(sprintf(
            "Possible corruption of gene names by excel: %s",
            paste(date_like, collapse = ", ")
        ))
    }
    #### Then fix with ####
    exp_CORRECTED <- exp
    yy <- HGNChelper::checkGeneSymbols(rownames(exp_CORRECTED),
        unmapped.as.na = TRUE
    )
    xx <- HGNChelper::checkGeneSymbols(rownames(exp_CORRECTED),
        unmapped.as.na = FALSE
    )
    numCorrected <- dim(xx[!xx$x == xx$Suggested.Symbol, ])[1]
    numBad <- sum(is.na(yy$Suggested.Symbol))
    messager(sprintf(
        "%s of %s gene symbols corrected.",
        numCorrected, dim(xx)[1]
    ), v = verbose)
    messager(sprintf(
        "%s of %s gene symbols cannot be mapped.",
        numBad, dim(xx)[1]
    ), v = verbose)
    newGnames <- xx$Suggested.Symbol
    exp_CORRECTED <- exp_CORRECTED[!duplicated(newGnames), ]
    newGnames <- newGnames[!duplicated(newGnames)]
    rownames(exp_CORRECTED) <- newGnames
    exp_CORRECTED <- to_sparse_matrix(
        exp = exp_CORRECTED,
        as_sparse = as_sparse,
        verbose = verbose
    )
    #### Now filter results in SE/SCE obj if inputted and return it ####
    if (SE_obj) {
        # Update all annotation datasets by replacing by corrected counts,
        # add in annotation and meta data if available
        # match the number of rows.
        SE_exp <- SE_exp[seq_len(nrow(exp_CORRECTED)), ]
        names(SummarizedExperiment::rowRanges(SE_exp)) <-
            rownames(exp_CORRECTED) # Update gene names
        SummarizedExperiment::assays(SE_exp) <- list(counts = exp_CORRECTED)
        return(SE_exp)
    }
    return(exp_CORRECTED)
}
