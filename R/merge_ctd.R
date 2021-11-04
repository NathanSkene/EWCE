#' Merge multiple CellTypeDataset references
#'
#' Import CellTypeDataset (CTD) references from a remote repository,
#' standardize each, and then merge into one CTD.
#'
#' Optionally, can return these as a merged \code{SingleCellExperiment}.
#'
#' @param CTD_list (Named) list of \code{CellTypeDatasets}.
#' @param save_dir The directory to save merged files in.
#' @param standardise_CTD Whether to run \code{standardise_ctd}.
#' @param as_SCE If \code{=T}, will return a the merged CTD
#' as a \code{SingleCellExperiment.}
#' Otherwise, will return a list of unmerged CTD.
#' @param gene_union Whether to take the gene union or intersection
#' when merging matrices (mean_exp,specificity, etc.).
#' @param merge_levels Which CTD levels you want to merge.
#' Can be a single value (e.g. \code{merge_levels=5})
#' or a list c(e.g. \code{merge_levels=c(1:5)}).
#' If some CTD don't have the same number of levels,
#' the maximum level depth available in that CTD will be used instead.
#' @param save_split_CTD Whether to save individual CTD files
#'  in the subdirectory \emph{standardized_CTD}.
#' @param save_split_SCE Whether to save individual SCE files
#'  in the subdirectory \emph{standardized_CTD_SCE}.
#' @param save_merged_SCE Save the final merged SCE object, or simply to return it.
#' @param force_new_quantiles If specificity quantiles matrix already exists, create new one.
#' @param numberOfBins Number of bins to compute specificity quantiles with.
#' @param as_sparse Convert matrices to sparse matrix.
#' @param as_DelayedArray Convert matrices to \code{DelayedArray}.
#' @param verbose Print messages.
#' @param ... Additional arguments to be passed to \code{standardise_ctd}.
#'
#' @return List of CellTypeDatasets or SingleCellExperiments.
#'
#' @examples
#' ## Let's pretend these are different CTD datasets
#' ctd1 <- ewceData::ctd()
#' ctd2 <- ctd1
#' CTD_list <- list(ctd1, ctd2)
#' SCE_merged <- EWCE::merge_ctd(
#'     CTD_list = CTD_list,
#'     as_SCE = TRUE,
#'     gene_union = TRUE
#' )
#' @export
#' @importFrom dplyr %>%
merge_ctd <- function(CTD_list,
                        save_dir = tempdir(),
                        standardise_CTD = FALSE,
                        as_SCE = FALSE,
                        gene_union = TRUE,
                        merge_levels = seq(1, 5),
                        save_split_SCE = FALSE,
                        save_split_CTD = FALSE,
                        save_merged_SCE = TRUE,
                        force_new_quantiles = FALSE,
                        numberOfBins = 40,
                        as_sparse = TRUE,
                        as_DelayedArray = FALSE,
                        verbose = TRUE,
                        ...) {
    #### Ensure it's actually a list ####
    if (all(is.null(names(CTD_list)))) {
        names(CTD_list) <-
            paste0("ctd", seq(1, length(CTD_list)))
    }
    #### Set up dirs ####
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    split_SCE_dir <- file.path(save_dir, "standardized_CTD_SCE")
    if (as_SCE) {
        dir.create(split_SCE_dir,
            showWarnings = FALSE, recursive = TRUE
        )
    }
    #### Standardise ####
    if (standardise_CTD) {
        CTD_list <- lapply(names(CTD_list), function(x) {
            message(x)
            standardise_ctd(
                ctd = CTD_list[[x]],
                ...
            )
        }) %>% `names<-`(names(CTD_list))
    }
    #### Split and save the standardized CTD files ####
    if (save_split_CTD) {
        split_CTD_dir <- file.path(save_dir, "standardized_CTD")
        dir.create(split_CTD_dir, showWarnings = TRUE, recursive = TRUE)
        CTD_paths <- lapply(names(CTD_list), function(x) {
            message(x)
            out_path <- file.path(split_CTD_dir, paste0(x, ".rds"))
            if (!file.exists(out_path)) {
                saveRDS(CTD_list[[x]], file = out_path)
            }
            return(out_path)
        })
    }
    #### Convert to nested list of SingleCellExperiments ####
    if (as_SCE) {
        messager("Converting CTD to merged SCEs", v = verbose)
        #### Convert to SCE ####
        SCE_lists <- lapply(names(CTD_list), function(x) {
            messager(x, v = verbose)
            sce <- ctd_to_sce(object = CTD_list[[x]])
            if (save_split_SCE) {
                saveRDS(sce, file.path(split_SCE_dir, paste0(x, ".rds")))
            }
            return(sce)
        }) %>% `names<-`(names(CTD_list))

        #### Split the standardized SCE files ####
        SCE_paths <- lapply(names(SCE_lists), function(x) {
            file.path(split_SCE_dir, paste0(x, ".rds"))
        }) %>% `names<-`(names(SCE_lists))
        #### Merge SCE ####
        SCE_merged <- merge_sce_list(
            SCE_lists = SCE_lists,
            parent_folder = split_SCE_dir,
            merge_levels = merge_levels,
            gene_union = gene_union,
            as_sparse = as_sparse,
            as_DelayedArray = as_DelayedArray,
            verbose = verbose
        )
        if (save_merged_SCE) {
            merged_dir <- file.path(save_dir, "merged")
            save_path <- paste0(
                merged_dir, "/CTD_SCE_merged.",
                if (gene_union) "union" else "intersect", ".rds"
            )
            dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)
            messager("Saving SCE_merged ==>", save_path, v = verbose)
            saveRDS(SCE_merged, save_path)
        }
        return(SCE_merged)
    } else {
        if (save_merged_SCE) {
            messager("+ Must set `as_SCE=TRUE` in order to merge CTD.",
                v = verbose
            )
            messager("+ Returning CTD_list.",
                v = verbose
            )
        }
        return(CTD_list)
    }
}
