# Use the rank norm transformation on specificity
rNorm <- function(ctdIN,
                  as_sparse = TRUE,
                  verbose = TRUE) {
    spec <- apply(ctdIN$specificity, 2, RNOmni::RankNorm)
    ctdIN$specificity <- to_sparse_matrix(
        exp = spec,
        as_sparse = as_sparse,
        verbose = verbose
    )
    return(ctdIN)
}
