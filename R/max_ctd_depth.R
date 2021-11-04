#' Get max CTD depth
#'
#' Get the maximum level depth from a list of CellTypeDataset objects.
#'
#' @param CTD_list A list of CellTypeDataset objects.
#'
#' @keywords internal
max_ctd_depth <- function(CTD_list) {
    max(unlist(lapply(CTD_list, length)))
}
