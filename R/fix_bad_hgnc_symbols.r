#' fix_bad_hgnc_symbols
#' - Given an expression matrix, wherein the rows are supposed to be HGNC 
#' symbols, find those symbols which are not official HGNC symbols, then
#' correct them if possible. Return the expression matrix with corrected 
#' symbols.
#' @param exp An expression matrix where the rows are MGI symbols
#' @param dropNonHGNC Boolean. Should symbols not recognised as HGNC symbols 
#' be dropped?
#' @return Returns the expression matrix with the rownames corrected and rows 
#' representing the same gene merged
#' @export
#' @import HGNChelper
#' @import biomaRt
#' @import ewceData
#' @import ExperimentHub
#' @importFrom AnnotationHub query 
fix_bad_hgnc_symbols <- function(exp, dropNonHGNC = FALSE) {

    # Check that exp is not some weird input format like those generated 
    # by readr functions
    if (!is(exp)[1] %in% c("matrix", "data.frame")) {
        stop("ERROR: exp must be either matrix or data.frame")
    }

    # First, find which gene symbols are not proper HGNC symbols
    #eh <- query(ExperimentHub::ExperimentHub(), "ewceData")
    #all_hgnc <- eh[["EH5371"]]
    all_hgnc <- ewceData::all_hgnc()
    not_HGNC <- rownames(exp)[!rownames(exp) %in% all_hgnc]
    print(sprintf("%s of %s are not proper HGNC symbols", 
                    length(unique(not_HGNC)), dim(exp)[1]))

    # If dropNonHGNC==TRUE then dropNonHGNC symbols
    if (dropNonHGNC == TRUE) {
        exp <- exp[!rownames(exp) %in% not_HGNC, ]
    }

    # First, check if any gene symbols have been corrupted by excel
    date_like <- not_HGNC[grep("SEP|MAR|FEB|DEC|Sep|Mar|Feb|Dec", not_HGNC)]
    if (length(date_like) > 0) {
        print(sprintf("Possible corruption of gene names by excel: %s", 
                        paste(date_like, collapse = ", ")))
        warning(sprintf("Possible corruption of gene names by excel: %s", 
                            paste(date_like, collapse = ", ")))
    }

    # Then fix with
    exp_CORRECTED <- exp
    yy <- checkGeneSymbols(rownames(exp_CORRECTED), unmapped.as.na = TRUE)
    xx <- checkGeneSymbols(rownames(exp_CORRECTED), unmapped.as.na = FALSE)
    numCorrected <- dim(xx[!xx$x == xx$Suggested.Symbol, ])[1]
    numBad <- sum(is.na(yy$Suggested.Symbol))
    print(sprintf("%s of %s gene symbols corrected", numCorrected, dim(xx)[1]))
    print(sprintf("%s of %s gene symbols cannot be mapped", numBad, dim(xx)[1]))
    newGnames <- xx$Suggested.Symbol
    exp_CORRECTED <- exp_CORRECTED[!duplicated(newGnames), ]
    newGnames <- newGnames[!duplicated(newGnames)]
    rownames(exp_CORRECTED) <- newGnames
    return(exp_CORRECTED)
}
