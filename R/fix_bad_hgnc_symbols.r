#' fix_bad_hgnc_symbols
#' - Given an expression matrix, wherein the rows are supposed to be HGNC
#' symbols, find those symbols which are not official HGNC symbols, then
#' correct them if possible. Return the expression matrix with corrected
#' symbols.
#' @param exp An expression matrix where the rows are HGNC symbols or an SingleCellExperiment (SCE) or other Ranged Summarized Experiment (SE) type object.
#' @param dropNonHGNC Boolean. Should symbols not recognised as HGNC symbols
#' be dropped?
#' @return Returns the expression matrix with the rownames corrected and rows
#' representing the same gene merged. If a SingleCellExperiment (SCE) or other
#' Ranged Summarized Experiment (SE) type object was inputted this will be
#' returned with the corrected expression matrix under counts.
#' @export
#' @import HGNChelper
#' @import biomaRt
#' @import ewceData
#' @import ExperimentHub
#' @importFrom AnnotationHub query
#' @importFrom SummarizedExperiment rowRanges assays
fix_bad_hgnc_symbols <- function(exp, dropNonHGNC = FALSE) {
    err_msg <- paste0("ERROR: 'exp' is null. It should be a numerical",
                      " matrix with the rownames being HGNC symbols.")
    # Check arguments
    if (is.null(exp)) {
        stop(err_msg)
    }
    #Check if input is an SCE or SE object
    if(methods::is(exp,"SummarizedExperiment")){
        #update exp to hold the counts from the SCE
        SE_exp<-exp
        if(!"counts" %in% names(SummarizedExperiment::assays(SE_exp)))
            stop(paste0("Please ensure counts is the assay name for your raw ",
                        "experiment data in your SE/SCE object"))
        exp <- SummarizedExperiment::assays(SE_exp)$counts
        #set boolean for later operations
        SE_obj<-TRUE
    } else {
        SE_obj<-FALSE
    }
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
    message(sprintf("%s of %s are not proper HGNC symbols",
                        length(unique(not_HGNC)), dim(exp)[1]))

    # If dropNonHGNC==TRUE then dropNonHGNC symbols
    if (dropNonHGNC == TRUE) {
        exp <- exp[!rownames(exp) %in% not_HGNC, ]
    }

    # First, check if any gene symbols have been corrupted by excel
    date_like <- not_HGNC[grep("SEP|MAR|FEB|DEC|Sep|Mar|Feb|Dec", not_HGNC)]
    if (length(date_like) > 0) {
        message(sprintf("Possible corruption of gene names by excel: %s",
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
    message(sprintf("%s of %s gene symbols corrected",
                        numCorrected, dim(xx)[1]))
    message(sprintf("%s of %s gene symbols cannot be mapped",
                        numBad, dim(xx)[1]))
    newGnames <- xx$Suggested.Symbol
    exp_CORRECTED <- exp_CORRECTED[!duplicated(newGnames), ]
    newGnames <- newGnames[!duplicated(newGnames)]
    rownames(exp_CORRECTED) <- newGnames
    #Now filter results in SE/SCE obj if inputted and return it
    if(SE_obj){
        #Update all annotation datasets by replacing by corrected counts, add
        # in annotation and meta data if available
        SE_exp <- SE_exp[seq_len(nrow(exp_CORRECTED)),] #match the number of rows
        names(SummarizedExperiment::rowRanges(SE_exp)) <-
            rownames(exp_CORRECTED) #Update gene names
        SummarizedExperiment::assays(SE_exp) <- list(counts=exp_CORRECTED)
        return(SE_exp)
    }
    return(exp_CORRECTED)
}
