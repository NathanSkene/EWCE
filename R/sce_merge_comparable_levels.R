# sce_merge_comparable_levels <- function(SCE_merged,
#                                         CTD_meta,
#                                         method="union"){
#   SCE_list <- lapply(unique(SCE_merged$level_1$batch), function(x){
#     print(x)
#     comp_lvl <- CTD_meta[CTD_meta$dataset==x, ]$comparison_level
#     print(comp_lvl)
#     sce <- SCE_merged[[comp_lvl]][,SCE_merged[[comp_lvl]]$batch==x]
#     for(ass in SummarizedExperiment::assayNames(sce)){
#       SummarizedExperiment::assay(sce, ass) <- as(  SummarizedExperiment::assay(sce, ass), "matrix")
#     }
#     if(ncol(sce)<30) warning(paste(x,": only contains",ncol(sce),"samples. This may cause problems with LIGER.\n\n"))
#     return(sce)
#   }) %>% `names<-`(unique(SCE_merged$level_1$batch))
#   sce_lvl <- sce_cbind(sce_list = SCE_list,
#                        method = method,
#                        batch_names = names(SCE_list),
#                        exprs = c("mean_exp","specificity","specificity_quantiles"))
#   return(sce_lvl)
# }
