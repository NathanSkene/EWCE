#' get.length.of.longest.string
#'
#' \code{get.length.of.longest.string} Given an array of strings (i.e. cellnames). Find the length of the longest
#' in terms of characters. This can be used to determine the width of PDF figures.
#'
#' @param x Array of names
#' @return Length of the longest name [integer]
#' @examples
#' get.length.of.longest.string(c("nathan","tom"))
#' @export
get.length.of.longest.string <- function(x){
    name_widths = strsplit(as.character(unique(x)),split=character(0))
    longest=0
    for(kk in 1:length(name_widths)){
        if(length(name_widths[[kk]])>longest){longest=length(name_widths[[kk]])}
    }
    return(longest)
}

