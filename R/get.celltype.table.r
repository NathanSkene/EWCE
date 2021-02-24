#' get.celltype.table
#'
#' \code{get.celltype.table} Generates a table that can be used for 
#' supplemenary tables of publications.
#' The table lists how many cells are associated with each celltype, the level 
#' of annotation, and the dataset from which it was generated.
#'
#' @param annot An annotation dataframe, which columns named 'level1class', 
#' 'level2class' and 'dataset_name'
#' @return A dataframe with columns 'name', 'level', 'freq' and 'dataset_name'
#' @examples
#' library(ewceData)
#' # See PrepLDSC.Rmd for origin of merged_ALLCELLS$annot
#' cortex_mrna <- cortex_mrna()
#' cortex_mrna$annot$dataset_name <- "cortex_mrna"
#' celltype_table <- get.celltype.table(cortex_mrna$annot)
#' @export
get.celltype.table <- function(annot) {
    err_msg <-paste0("Error: annotation dataframe must have columns for",
                        " 'level1class', 'level2class' and 'dataset_name")
    if (is.null(annot$level1class) | is.null(annot$level1class) | 
            is.null(annot$dataset_name)) {
        stop(err_msg)
    }
    annot$cell_id <- as.character(annot$cell_id)
    annot$level1class <- as.character(annot$level1class)
    annot$level2class <- as.character(annot$level2class)
    annot$dataset_name <- as.character(annot$dataset_name)

    compress_name <- function(x) {
        b <- gsub(" ", "", x)
        c <- gsub("/", "and", b)
        return(c)
    }
    l1_datasets <- unique(annot[, c("level1class", "dataset_name")])
    colnames(l1_datasets)[1] <- "name"
    l1_datasets[, 1] <- compress_name(l1_datasets[, 1])
    l2_datasets <- unique(annot[, c("level2class", "dataset_name")])
    colnames(l2_datasets)[1] <- "name"
    l2_datasets[, 1] <- compress_name(l2_datasets[, 1])

    # Find all L2 celltypes which contribute to each L1 celltype
    l1_contribs <- rep("", length(l1_datasets$name))
    names(l1_contribs) <- l1_datasets$name
    for (l1 in l1_datasets$name) {
        l1_contribs[l1] <- 
            paste(unique(annot[
                compress_name(annot$level1class) == l1, ]$level2class), 
                    collapse = ", ")
    }
    l1_contribs <- data.frame(name = names(l1_contribs), 
                                contributing_l2_celltypes = l1_contribs, 
                                stringsAsFactors = FALSE)
    l1_datasets$contributing_l2_celltypes <- 
        l1_contribs$contributing_l2_celltypes
    l2_datasets$contributing_l2_celltypes <- "NA"

    # Calculate frequency of level 1 celltypes
    freqs <- rep(0, dim(l1_datasets)[1])
    for (r in seq_len(dim(l1_datasets)[1])) {
        freqs[r] <- 
            dim(annot[
                compress_name(annot$level1class) == l1_datasets[r, ]$name & 
                    annot$dataset_name == l1_datasets[r, ]$dataset_name, ])[1]
    }
    l1_datasets$number_of_cells <- freqs

    # Calculate frequency of level 2 celltypes
    freqs <- rep(0, dim(l2_datasets)[1])
    for (r in seq_len(dim(l2_datasets)[1])) {
        freqs[r] <- 
            dim(annot[
                compress_name(annot$level2class) == l2_datasets[r, ]$name & 
                    annot$dataset_name == l2_datasets[r, ]$dataset_name, ])[1]
    }
    l2_datasets$number_of_cells <- freqs

    l1_datasets$annotation_level <- 1
    l2_datasets$annotation_level <- 2
    l1_ct <- l1_datasets[, c("name", "annotation_level", "dataset_name", 
                                "number_of_cells", "contributing_l2_celltypes")]
    l2_ct <- l2_datasets[, c("name", "annotation_level", "dataset_name", 
                                "number_of_cells", "contributing_l2_celltypes")]
    celltype_table <- rbind(l1_ct, l2_ct)

    return(celltype_table)
}
