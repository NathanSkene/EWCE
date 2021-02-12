#' drop.uninformative.genes
#'
#' \code{drop.uninformative.genes} drops genes from an SCT expression matrix 
#' if they do not significantly vary between any celltypes.
#' Makes this decision based on use of an ANOVA (implemented with Limma). If 
#' the F-statistic for variation amongst type2 annotations
#' is less than a strict p-threshold, then the gene is dropped.
#'
#' @param exp Expression matrix with gene names as rownames.
#' @param level2annot Array of cell types, with each sequentially corresponding 
#' a column in the expression matrix
#' @return exp Expression matrix with gene names as rownames.
#' @examples
#' data("cortex_mrna")
#  # Use only a subset of genes to keep the example quick
#' cortex_mrna$exp <- cortex_mrna$exp[1:300, ] 
#' exp2 <- drop.uninformative.genes(exp = cortex_mrna$exp, 
#'     level2annot = cortex_mrna$annot$level2class)
#' @export
#' @import limma
#' @import stats
#' @importFrom methods is
drop.uninformative.genes <- function(exp, level2annot) {
    if (is(exp[1, 1])[1] == "character") {
        exp <- as.matrix(exp)
        storage.mode(exp) <- "numeric"
    }
    level2annot <- as.character(level2annot)
    summed <- apply(exp, 1, sum)
    exp <- exp[summed != 0, ]
    mod <- model.matrix(~level2annot)
    fit <- lmFit(exp, mod)
    eb <- eBayes(fit)
    pF <- p.adjust(eb$F.p.value, method = "BH")
    exp <- exp[pF < 0.00001, ]
    return(exp)
}
