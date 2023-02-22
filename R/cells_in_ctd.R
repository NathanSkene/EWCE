cells_in_ctd <- function(ctdIN, 
                         cells) {
    if (sum(!cells %in% colnames(ctdIN$specificity) == 0)) {
        return(1)
    } else {
        return(0)
    }
} 