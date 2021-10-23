merge_sce_list <- function(SCE_lists = NULL,
                           parent_folder = NULL,
                           pattern = ".rds",
                           merge_levels = seq(1, 5),
                           gene_union = TRUE,
                           as_sparse = TRUE,
                           as_DelayedArray = TRUE,
                           verbose = TRUE) {
    #### If SCE_lists hasn't been created yet, create it ###
    if (is.null(SCE_lists) && (!is.null(parent_folder))) {
        sce_files <- list.files(
            parent_folder,
            pattern = pattern,
            ignore.case = TRUE,
            full.names = TRUE
        )
        messager(paste(
            length(sce_files),
            "SingleCellExperiment files found."
        ),
        v = verbose
        )
        #### Import files ####
        SCE_lists <- lapply(sce_files, function(x) {
            dataset <- gsub(pattern, "", basename(x))
            print(dataset)
            sce <- readRDS(x)
            return(sce)
        }) %>% `names<-`(gsub(pattern, "", basename(sce_files)))
    }
    #### Merge SCE_lists into one SCE ####
    max_depth <- max_ctd_depth(SCE_lists)
    merge_levels <- seq(
        min(merge_levels),
        min(max(merge_levels), max_depth)
    )

    SCE_merged <- lapply(merge_levels, function(lvl) {
        messager("Merging SCE at level:", lvl, v = verbose)
        sce.lvl <- sce_lists_apply(SCE_lists,
            level = lvl,
            as_matrix = TRUE
        )
        sce_cbind(
            sce_list = sce.lvl,
            method = if (gene_union) "union" else "intersect",
            batch_names = names(sce.lvl),
            exprs = c("mean_exp", "specificity", "specificity_quantiles")
        )
    }) %>% `names<-`(paste("level", merge_levels, sep = "_"))
    #### Convert back to sparseDM ####
    SCE_merged <- sce_merged_apply(SCE_merged,
        as_sparse = as_sparse,
        as_DelayedArray = as_DelayedArray
    )
    return(SCE_merged)
}
