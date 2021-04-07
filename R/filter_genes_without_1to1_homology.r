#' filter_genes_without_1to1_homolog
#'
#' \code{filter_genes_without_1to1_homolog} Takes the filenames of
#' celltype_data files, loads them, and drops any genes
#' which don't have a 1:1 homolog based on biomart. The new files are saved to
#' disc, appending '_1to1only' to the file
#' tag of the previous file.
#'
#' @param filenames Array of filenames for sct_data saved as .Rda files
#' @return Array of filenames included the ones with only 1:1 homologs
#' @examples
#' library(ewceData)
#' # Load the single cell data
#' cortex_mrna <- cortex_mrna()
#' expData <- cortex_mrna$exp
#' expData <- expData[1:100, ] # Use only a subset to keep the example quick
#' l1 <- cortex_mrna$annot$level1class
#' l2 <- cortex_mrna$annot$level2class
#' annotLevels <- list(level1class = l1, level2class = l2)
#' fNames_ALLCELLS <- generate_celltype_data(exp = expData,
#'     annotLevels = annotLevels,
#'     groupName = "allKImouse")
#' fNames_ALLCELLS <- filter_genes_without_1to1_homolog(fNames_ALLCELLS)
#' @export
#' @import utils
#' @import stringr
#' @import ewceData
#' @import ExperimentHub
#' @importFrom AnnotationHub query
filter_genes_without_1to1_homolog <- function(filenames) {
    newFilenames <- filenames
    mouse_to_human_homologs <- ewceData::mouse_to_human_homologs()
    #quicker to load with ID
    #eh <- query(ExperimentHub::ExperimentHub(), "ewceData")
    #mouse_to_human_homologs <- eh[["EH5367"]]
    orthologs <- mouse_to_human_homologs
    mgi_1to1 <- orthologs$MGI.symbol
    hgnc_1to1 <- orthologs$HGNC.symbol
    for (ff in filenames) {
        load(ff)
        sct_genes <- rownames(ctd[[1]]$mean_exp)
        # If it's a m
        if (sum(sct_genes %in% orthologs$MGI.symbol) >
                sum(sct_genes %in% orthologs$HGNC.symbol)) {
            symbol_1to1_in_sct <- mgi_1to1[mgi_1to1 %in% sct_genes]
        } else {
            symbol_1to1_in_sct <- hgnc_1to1[hgnc_1to1 %in% sct_genes]
        }
        ctd[[1]]$mean_exp <- ctd[[1]]$mean_exp[symbol_1to1_in_sct, ]
        ctd[[1]]$specificity <- ctd[[1]]$specificity[symbol_1to1_in_sct, ]
        # ff2 = gsub("_level","_1to1only_level",ff)
        ff2 <- gsub("\\.rda", "_1to1only\\.rda", ff)
        save(ctd, file = ff2)
        newFilenames <- c(newFilenames, ff2)
    }
    return(newFilenames)
}
