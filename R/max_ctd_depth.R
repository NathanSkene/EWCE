max_ctd_depth <- function(CTD_list) {
    max(unlist(lapply(CTD_list, length)))
}
