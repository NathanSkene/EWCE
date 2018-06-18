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


convert_new_ewce_to_old <- function(ctd,lvl){
    celltype_data=list()
    celltype_data[[1]]=list()
    celltype_data[[1]]$cell_dists  = ctd[[lvl]]$specificity
    celltype_data[[1]]$all_scts    = ctd[[lvl]]$mean_exp
    #celltype_data[[2]]=list()
    #celltype_data[[2]]$cell_dists = ctd[[2]]$specificity
    #celltype_data[[2]]$all_scts   = ctd[[2]]$mean_exp
    return(celltype_data)
}
