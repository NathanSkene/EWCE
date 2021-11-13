#' Fix celltype names
#'
#' Make sure celltypes don't contain characters that could interfere with
#' downstream analyses. For example, the R package
#' \href{https://github.com/neurogenomics/MAGMA_Celltyping}{MAGMA.Celltyping}
#' cannot have spaces in celltype names because spaces are used as a delimiter
#' in later steps.
#'
#' @param celltypes Character vector of celltype names.
#' @param replace_chars Regex string of characters to replace
#'  with "_" when renaming columns.
#'
#' @returns Fixed celltype names.
#'
#' @export
#' @examples
#' ct <- c("microglia", "astryocytes", "Pyramidal SS")
#' ct_fixed <- fix_celltype_names(celltypes = ct)
fix_celltype_names <- function(celltypes,
                               replace_chars = "[-]|[.]|[ ]|[//]|[\\/]") {
    if (is.null(celltypes)) {
        return(NULL)
    }
    celltypes <- gsub(replace_chars, "_", celltypes)
    # Remove repeat __
    celltypes <- gsub("[_+]", "_", celltypes)
    return(celltypes)
}
