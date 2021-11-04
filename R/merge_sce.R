#' Merge multiple \code{SingleCellExperiment} objects
#'
#' Merge several \code{SingleCellExperiment} (SCE) objects from
#' different batches/experiments.
#' Extracted from the
#' \href{https://bioconductor.org/packages/release/bioc/html/scMerge.html}{
#' scMerge} package.
#'
#' @param sce_list A list contains the \code{SingleCellExperiment}
#' Object from each batch.
#' @param method A string indicates the method of combining the
#' gene expression matrix, either \code{union} or \code{intersect}.
#'  Default to \code{intersect}. \code{union} only supports matrix class.
#' @param cut_off_batch A numeric vector indicating the cut-off for
#'  the proportion of a gene is expressed within each batch.
#' @param cut_off_overall A numeric vector  indicating the cut-off for
#' the proportion of a gene is expressed overall data.
#' @param use_assays A string vector indicating the expression matrices
#' to be combined.
#' The first assay named will be used to determine the proportion of zeros.
#' @param colData_names A string vector indicating the \code{colData}
#'  that are combined.
#' @param batch_names A string vector indicating the batch names for
#'  the output SCE object.
#' @param verbose Print messages.
#'
#' @return A \code{SingleCellExperiment} object with the list of SCE
#' objects combined.
#'
#' @source
#' \href{https://bioconductor.org/packages/release/bioc/html/scMerge.html}{
#' scMerge}.
#' @author Yingxin Lin (modified by Brian Schilder)
#'
#' @importFrom DelayedArray rowMeans cbind rbind
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' 
#' @export
#' @examples
#' ctd <- ewceData::ctd()
#' sce_list <- EWCE::ctd_to_sce(object = ctd)
#' sce_combine <- merge_sce(sce_list = sce_list) 
merge_sce <- function(sce_list,
                      method = "intersect",
                      cut_off_batch = 0.01,
                      cut_off_overall = 0.01,
                      use_assays = NULL,
                      colData_names = NULL,
                      batch_names = NULL,
                      verbose = TRUE) {
    #### Check args ####
    if(!method[1] %in% c("intersect","union")){
        stop("method must be one of: 'intersect', 'union'")
    }
    #### Get shared assays ####
    all_assays <- lapply(sce_list, function(sce) {
        names(SummarizedExperiment::assays(sce))
    } )
    shared_assays <- Reduce(intersect,all_assays)
    use_assays <- use_assays 
    if(is.null(use_assays) || (!use_assays %in% shared_assays)){
        use_assays <- shared_assays
    } 
    messager(
        "The assay", paste0("'",use_assays[1],"'"),
        "will be used to determine the",
        "proportion of zeroes for each batch.",
        v=verbose
    ) 
    
    if (method == "union" &&
        !is.matrix(SummarizedExperiment::assay(sce_list[[1]],
                                               use_assays[1]))) {
        stop_msg <- paste(
            "The union method only supports",
            "matrix class in the sce object \n "
        )
        stop(stop_msg)
    }
    n_batch <- length(sce_list)
    zero_list <- lapply(
        sce_list,
        function(x) {
            DelayedArray::rowMeans(
                SummarizedExperiment::assay(x, use_assays[1]) == 0
            )
        }
    )

    expressed_list <- lapply(zero_list, function(x) {
        names(which(x <=
            (1 - cut_off_batch)))
    })
    for (i in seq_len(n_batch)) {
        sce_list[[i]] <- sce_list[[i]][expressed_list[[i]], ]
    }
    
    #### Take gene intersect ####
    if (method == "intersect") {
        keep <- Reduce(intersect, expressed_list)
        sce_list <- lapply(sce_list, function(x) x[keep, ])
        assay_list <- list()
        for (i in seq_len(length(use_assays))) {
            assay_list[[i]] <- do.call(cbind, lapply(
                sce_list,
                function(y) SummarizedExperiment::assay(y, use_assays[i])
            ))
        }
        names(assay_list) <- use_assays
        colData_list <- do.call(
            DelayedArray::rbind,
            lapply(sce_list, function(y) {
                SummarizedExperiment::colData(y)[, colData_names, drop = FALSE]
            })
        )
        sce_combine <- SingleCellExperiment::SingleCellExperiment(
            assay = assay_list,
            colData = colData_list
        )
    }
    #### Take gene union ####
    if (method == "union") {
        keep <- Reduce(union, expressed_list)
        assay_list <- list()
        for (i in seq_len(length(use_assays))) {
            assay_list[[i]] <- do.call(
                DelayedArray::cbind,
                lapply(sce_list, function(x) {
                    mat <- matrix(0,
                        nrow = length(keep), ncol = ncol(x),
                        dimnames = list(keep, colnames(x))
                    )
                    mat[rownames(x), ] <-
                        SummarizedExperiment::assay(x, use_assays[i])
                    return(mat)
                })
            )
        }
        names(assay_list) <- use_assays
        colData_list <- do.call(rbind, lapply(sce_list, function(y) {
            SummarizedExperiment::colData(y)[,
                colData_names,
                drop = FALSE
            ]
        }))

        sce_combine <- SingleCellExperiment::SingleCellExperiment(
            assay = assay_list,
            colData = colData_list
        ) 
        zero_cbind <- DelayedArray::rowMeans(
            SummarizedExperiment::assay(sce_combine, use_assays[1]
                                        ) == 0
        )
        sce_combine <- sce_combine[names(zero_cbind[zero_cbind <=
            (1 - cut_off_overall)]), ]
    }
    #### Name batches ####
    if (is.null(batch_names)) {
        batch_names <- paste0("batch", seq_len(n_batch))
    }
    sce_combine$batch <- rep(batch_names, unlist(lapply(
        sce_list,
        ncol
    )))
    return(sce_combine)
}
