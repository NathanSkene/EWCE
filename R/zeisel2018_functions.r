#
#
# search_zeisel2018_celltypes <- function(target_celltype="purkinje",
#                                         verbose=T){
#     z18.celltypes <- ewceData::zeisel2018.celltypes()
#     z18.terms <- z18.celltypes[grep(target_celltype,z18.celltypes$Description,ignore.case = T),]$`Cluster name` |> unique()
#     if(verbose){
#         print(paste(length(z18.terms),"matching celltype(s)",
#                     "identified in zeisel2018."))
#     }
#     return(z18.terms)
# }
#
# translate_zeisel2018 <- function(sig_results){
#     print("+ Translating CellTypes from zeisel2018.")
#     z18.celltypes <- ewceData::zeisel2018.celltypes()
#
#     sig_results$CellType <- gsub("mouse.zeisel2018.","",sig_results$CellType)
#     sig_results2 <- merge(sig_results,
#                           z18.celltypes[,c("Cluster name","Description")],
#                           by.x="CellType", by.y="Cluster name",
#                           all.x = T,
#                           all.y = F)
#     # Make unique IDs
#     sig_results2$ClusterName <- sig_results2$CellType
#     sig_results2$CellType <- paste(gsub(" ","_",gsub(", ",".",sig_results2$Description)),
#                                    sig_results2$ClusterName,sep = ".")
#
#     row.names(sig_results2) <- make.unique(sig_results2$CellType)
#     return(sig_results2)
# }
#
#
# multilevel_zeisel2018_search <- function(ctd.zeisel2018,
#                                          target_celltype){
#     z18.terms <- lapply(1:length(ctd.zeisel2018), function(lvl){
#         print(lvl)
#         grep(target_celltype,
#              colnames(ctd.zeisel2018[[lvl]]$mean_exp), ignore.case = T, value = T)
#         # return(unique(colnames(ctd.zeisel2018[[lvl]]$mean_exp)))
#     }) |> `names<-`(names(ctd.zeisel2018))
# }
#
#
# get_zeisel2018_markers <- function(ctd.zeisel2018,
#                                    level=5,
#                                    target_celltype="purkinje",
#                                    use_specificity=T,
#                                    use_expression=T,
#                                    max_specificity_genes=NULL,
#                                    max_expression_genes=NULL,
#                                    translate_target_celltype=T){
#     specQ <- ctd.zeisel2018[[level]]$specificity_quantiles
#     spec <- ctd.zeisel2018[[level]]$specificity
#     mean_exp <- ctd.zeisel2018[[level]]$mean_exp
#     if(translate_target_celltype){
#         z18.terms <- search_zeisel2018_celltypes(target_celltype=target_celltype)
#         target_cols <- grep(paste(z18.terms,collapse = "|"),colnames(specQ))
#     }else {
#         # Get an exact match (instead of substring) because the level 5 labels are too imprecise
#         z18.terms <- target_celltype
#         target_cols <- colnames(mean_exp)[colnames(mean_exp) == paste0("mouse.zeisel2018.",target_celltype)]
#     }
#     if(!any(length(target_cols))) stop("No matching columns identified.")
#
#
#     #### Create mean expression quantiles matrix ####
#     mean_exp_q <- sapply(colnames(mean_exp), function(x){
#         bin_columns_into_quantiles(mean_exp[,x], 41)
#     }) |> `row.names<-`(row.names(mean_exp)) |>
#         `colnames<-`(colnames(mean_exp)) |>
#         as("sparseMatrix")
#     ## Sometimes you have to limit the number of genes or else everything is enriched
#     celltype_specific <- (sapply(target_cols, function(x){
#         g1 <- row.names(specQ)[specQ[,x]==max(specQ[,x], na.rm = T)]
#         if(!is.null(max_specificity_genes)){
#             top_genes <- data.frame(specificity=spec[,x], gene = row.names(spec) ) |>
#                 dplyr::arrange(desc(specificity))
#             g1 <- top_genes$gene[1:min(length(g1),max_specificity_genes, na.rm = T)]
#         }
#         return(g1)
#     }) |> reshape2::melt())$value |> unique()
#     print(paste("+",length(celltype_specific),"celltype_specific markers identified."))
#
#     celltype_expressed <- (sapply(target_cols, function(x){
#         g2 <- row.names(mean_exp_q)[ mean_exp_q[,x]==max(mean_exp_q[,x],na.rm = T)]
#         if(!is.null(max_expression_genes)) g2 <- g2[1:max_expression_genes]
#         return(g2)
#     }) |> reshape2::melt())$value |> unique()
#     print(paste("+",length(celltype_expressed),"celltype_expressed markers identified."))
#
#     ### Return different genes sets depending on arguments ####
#     if(use_specificity & use_expression==F){
#         print("+ Returning celltype_specific markers only.")
#         return(celltype_specific)
#     }
#     if(use_expression & use_specificity==F){
#         print("+ Returning celltype_expressed markers only.")
#         return(celltype_expressed)
#     }
#     if(use_specificity & use_expression){
#         celltype_markers <- intersect(celltype_specific, celltype_expressed)
#         print(paste("+ Returning the intersection of celltype_specific & celltype_expressed markers:",length(celltype_markers),"genes."))
#         return(celltype_markers)
#     }
#     if(use_specificity==F & use_expression==F){
#         stop("At least one of these must be set to TRUE:",
#              "use_specificity, use_expression")
#     }
# }
