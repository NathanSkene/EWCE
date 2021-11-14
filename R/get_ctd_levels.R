#' Get the names of CellTypeDataset levels
#'
#' Returns the level names of a CellTypeDataset. If none are available,
#' will instead return a vector of numbers (one number per level).
#'
#' @param ctd CellTypeDataset.
#' @param max_only Only return the level with the greatest depth
#'  (e.g. \code{"level3"} in \code{c("level1","level2","level3")}).
#'  
#' @return List of levels in \code{ctd}.
#'  
#' @keywords internal
get_ctd_levels <- function(ctd,
                           max_only = FALSE) {
    # This is necessary in case further meta-data such as $name is used
    if (!is.null(names(ctd))) {
        lvls <- names(ctd)
    } else {
        lvls <- seq(1, length(ctd))
    }
    if (max_only) {
        return(max(lvls))
    } else {
        return(lvls)
    }
}
