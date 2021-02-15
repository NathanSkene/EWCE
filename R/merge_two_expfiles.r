#' merge_two_expfiles
#'
#' \code{merge_two_expfiles} Used to combine two single cell type datasets
#'
#' @param exp1 Numerical expression matrix for dataset1 with row for each gene 
#' and column for each cell. Row names are MGI/HGNC gene symbols. Column names 
#' are cell IDs which can be cross referenced against the annot data frame.
#' @param exp2 Numerical expression matrix for dataset2 with row for each gene 
#' and column for each cell. Row names are MGI/HGNC gene symbols. Column names 
#' are cell IDs which can be cross referenced against the annot data frame.
#' @param annot1 Annotation data frame for dataset1 which contains three 
#' columns at least: cell_id, level1class and level2class
#' @param annot2 Annotation data frame for dataset2 which contains three 
#' columns at least: cell_id, level1class and level2class
#' @param name1 Name used to refer to dataset 1. Leave blank if it's already a 
#' merged dataset.
#' @param name2 Name used to refer to dataset 2. Leave blank if it's already a 
#' merged dataset.
#' @return List containing merged exp and annot
#' @examples
#' # See the EWCE vignette for further explanation
#' # Download the hypothalamus data and unzip
#' if (!file.exists("GSE74672_expressed_mols_with_classes.xlsx")) {
#'     file_link1 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74672/"
#'     file_link2 <- "suppl/GSE74672_expressed_mols_with_classes.xlsx.gz"
#'     download.file(paste0(file_link1, file_link2),
#'                    destfile = "GSE74672_expressed_mols_with_classes.xlsx.gz"
#'     )
#'     system("gunzip GSE74672_expressed_mols_with_classes.xlsx.gz")
#' }
#' # Read in the hypothalamus data -
#' # first 100 genes and subset of samples only to speed up computation
#' hypo_dat <- readxl::read_excel("GSE74672_expressed_mols_with_classes.xlsx",
#'     range="A1:CBD113")
#' # Extract the expression data, gene symbols and annotation data
#' exp <- data.matrix(hypo_dat[12:dim(hypo_dat)[1], 2:dim(hypo_dat)[2]])
#' rownames(exp) <- data.frame(hypo_dat[12:dim(hypo_dat)[1], 1])[, 1]
#' level1class <- data.frame(level1class = t(hypo_dat[1, 2:dim(hypo_dat)[2]]), 
#'                               stringsAsFactors = FALSE)[, 1]
#' level2class <- data.frame(leve2class = t(hypo_dat[2, 2:dim(hypo_dat)[2]]), 
#'                               stringsAsFactors = FALSE)[, 1]
#' cell_id <- colnames(hypo_dat)[2:dim(hypo_dat)[2]]
#' hypo_annot <- data.frame(
#'     cell_id = cell_id, level1class = level1class,
#'     level2class = level2class, stringsAsFactors = FALSE
#' )
#' # Drop the glia and unclassified cells(which don't have level 2  annotations)
#' hypo_annot <- 
#'     hypo_annot[!is.na(hypo_annot$level2class) & 
#'         !hypo_annot$level2class == "uc", ]
#' hypo_exp <- exp[, hypo_annot$cell_id]
#' # Make the celltype names more aesthetically pleasing
#' hypo_annot$level2class <- gsub(",", ";", hypo_annot$level2class)
#' hypo_annot$level1class[grep("Oxt;|^Avp", hypo_annot$level2class)] <-
#'     "Oxytocin / Vasopressin Expressing Neurons"
#' hypo_annot$level1class[grep("^Th;|^Dopamine", hypo_annot$level2class)] <-
#'     "Hypothalamic Dopaminergic Neurons"
#' hypo_annot$level1class[grepl("^Vglut2|^Trh|^Qrfp|^Hcrt|^Pmch|^Adcyap1|
#'                                    ^Npvf|^Ghrh|^Hmit|
#'                                    ^Nms|^Vip;|^Per2|Tnr$|^Gad-low;Gnrh",
#'     hypo_annot$level2class
#' ) & grepl("neurons", hypo_annot$level1class)] <-
#'     "Hypothalamic Glutamatergic Neurons"
#' hypo_annot$level1class[grepl(
#'     "GABA|^Sst|^Crh|^Npy|^Pomc|^Galanin|^Otof|Pnoc$|^Calcr-high",
#'     hypo_annot$level2class
#' ) & grepl("^neurons$", hypo_annot$level1class)] <-
#'     "Hypothalamic GABAergic Neurons"
#' hypo_annot$level2class[hypo_annot$level2class != ""] <-
#'     sprintf("Hypothalamic %s Neuron", 
#'         hypo_annot$level2class[hypo_annot$level2class != ""])
#' # Fix bad MGI symbols
#' hypo_exp_CORRECTED <- fix.bad.mgi.symbols(hypo_exp)
#' # Merge hypothalamus data with the cortex dataset
#' merged_KI <- merge_two_expfiles(
#'     exp1 = hypo_exp_CORRECTED, exp2 = cortex_mrna$exp,
#'     annot1 = hypo_annot, annot2 = cortex_mrna$annot,
#'     name1 = "Hypothalamus (KI)", name2 = "Cortex/Hippo (KI)"
#' )
#' @export
merge_two_expfiles <- function(exp1, exp2, annot1, annot2, 
                                name1 = "", name2 = "") {
    # genes_intersect   = intersect(rownames(exp1),rownames(exp2))
    # genes_exp1_unique = setdiff(rownames(exp1),rownames(exp2))
    # genes_exp2_unique = setdiff(rownames(exp2),rownames(exp1))
    all_genes <- unique(c(rownames(exp2), rownames(exp1)))

    # Check for correct column headers in annot data tables
    if (sum(c("cellid") %in% colnames(annot1)) == 1) {
        colnames(annot1)[colnames(annot1) == "cellid"] <- "cell_id"
    }
    if (sum(c("cellid") %in% colnames(annot2)) == 1) {
        colnames(annot2)[colnames(annot2) == "cellid"] <- "cell_id"
    }
    err_msg <- paste0("ERROR: annot1 doesn't have either cell_id,",
                        " level1class or level2class columns")
    err_msg2 <- paste0("ERROR: annot1 doesn't have either cell_id,",
                        " level1class or level2class columns")
    if (!sum(c("cell_id", "level1class", "level2class") %in% 
                colnames(annot1)) == 3) {
        stop(err_msg)
    }
    if (!sum(c("cell_id", "level1class", "level2class") %in% 
                colnames(annot2)) == 3) {
        stop(err_msg2)
    }

    # If one of the exp matrices is really a matrix (not a dataframe) 
    # the code won't work... so force conversion
    if (is(exp1)[1] == "matrix") {
        nn <- colnames(exp1)
        rr <- rownames(exp1)
        exp1 <- data.frame(exp1, stringsAsFactors = FALSE)
        colnames(exp1) <- nn
        rownames(exp1) <- rr
    }
    if (is(exp2)[1] == "matrix") {
        nn <- colnames(exp2)
        rr <- rownames(exp2)
        exp2 <- data.frame(exp2, stringsAsFactors = FALSE)
        colnames(exp2) <- nn
        rownames(exp2) <- rr
    }


    # Merge the expression matrices, setting undetected genes to 0
    exp1b <- exp1[all_genes, ]
    exp1b[is.na(exp1b)] <- 0
    rownames(exp1b) <- all_genes
    exp2b <- exp2[all_genes, ]
    exp2b[is.na(exp2b)] <- 0
    rownames(exp2b) <- all_genes
    exp <- cbind(exp1b, exp2b)

    # Ensure important annotation columns are not stored as factors
    annot1$cell_id <- as.character(annot1$cell_id)
    annot2$cell_id <- as.character(annot2$cell_id)
    annot1$level1class <- as.character(annot1$level1class)
    annot1$level2class <- as.character(annot1$level2class)
    if (is.null(annot1$dataset_name)) {
        annot1$dataset_name <- name1
    }
    annot2$level1class <- as.character(annot2$level1class)
    annot2$level2class <- as.character(annot2$level2class)
    if (is.null(annot2$dataset_name)) {
        annot2$dataset_name <- name2
    }
    keepTISSUE <- FALSE
    if (("tissue" %in% colnames(annot1)) & ("tissue" %in% colnames(annot2))) {
        keepTISSUE <- TRUE
        annot1$tissue <- as.character(annot1$tissue)
        annot2$tissue <- as.character(annot2$tissue)
    }

    # Setup new annotation dataframe
    cell_id <- c(annot1$cell_id, annot2$cell_id)
    level1class <- c(annot1$level1class, annot2$level1class)
    level2class <- c(annot1$level2class, annot2$level2class)
    dataset_name <- c(as.character(annot1$dataset_name), 
                        as.character(annot2$dataset_name))
    if (!keepTISSUE) {
        annot <- data.frame(cell_id = cell_id, level1class = level1class, 
                                level2class = level2class, 
                                dataset_name = dataset_name)
    } else {
        tissue <- c(annot1$tissue, annot2$tissue)
        annot <- data.frame(cell_id = cell_id, level1class = level1class, 
                                level2class = level2class, 
                                dataset_name = dataset_name, tissue = tissue)
    }


    # Drop expression data without annotation data
    numMissingAnnot <- dim(exp)[2] - sum(annot$cell_id %in% colnames(exp))
    if (numMissingAnnot > 0) {
        txt <- paste0("Warning: %s cells are missing annotation data",
                        " and have been dropped")
        print(sprintf(txt, numMissingAnnot))
        exp <- exp[, as.character(annot$cell_id)]
    }

    return(list(exp = exp, annot = annot))
}
