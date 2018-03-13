convert_old_ewce_to_new <- function(level1,level2){
    load(level1)
    ctd=list()
    ctd[[1]]=list()
    ctd[[1]]$specificity = celltype_data[[1]]$cell_dists
    ctd[[1]]$mean_exp    = celltype_data[[1]]$all_scts
    load(level2)
    ctd[[2]]=list()
    ctd[[2]]$specificity = celltype_data[[1]]$cell_dists
    ctd[[2]]$mean_exp    = celltype_data[[1]]$all_scts
    return(ctd)
}