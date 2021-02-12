#' load.linnarsson.sct.data
#'
#' \code{load.linnarsson.sct.data} Function which was specifically written for 
#' extracting the expression and annotation data from the publically 
#' downloadable Linnarsson lab Cortex/Hippocampus dataset
#'
#' @param fName Path to the downloaded file
#' @return A list containing 'exp' and 'annot'.
#' \itemize{
#'   \item \code{exp}: Numerical matrix with row for each gene and column for 
#'   each cell. Row names are MGI/HGNC gene symbols. Column names are cell IDs 
#'   which can be cross referenced against the annot data frame.
#'   \item \code{annot}: Annotation data frame which must contain at least 
#'   three columns: cell_id, level1class and level2class. Can also contain 
#'   other columns.
#' }
#' @examples
#' # Load the single cell data
#' linURL <- "goo.gl/r5Y24y"
#' download.file(linURL, destfile = "expression_mRNA_17-Aug-2014.txt")
#' path <- "expression_mRNA_17-Aug-2014.txt"
#' cortex_mrna <- load.linnarsson.sct.data(path)
#' @export
#' @import utils
load.linnarsson.sct.data <- function(fName) {
    all_dat <- read.csv(fName, sep = "\t", stringsAsFactors = FALSE)
    # Check that data is formatted as expected
    corr_fmt<-paste0("group #,total mRNA mol,well,sex,age,",
                        "diameter,cell_id,level1class,level2class")
    if (paste(all_dat[seq_len(9), 2], collapse = ",") != corr_fmt) {
        stop("ERROR: annotation rows are not written as expected")
    }
    # A blank row seperates the annotation from expression data...find that row
    expDataStarts <- 1
    for (i in seq_len(30))
    {
        if (length(unique(unlist(all_dat[i, ]))) == 1 & 
                unique(unlist(all_dat[i, ]))[1] == "") {
            expDataStarts <- i + 1
        }
    }
    # Get the expression data and convert to numerical matrix
    zeisel_exp <- all_dat[expDataStarts:dim(all_dat)[1], 3:dim(all_dat)[2]]
    zeisel_exp <- as.matrix(zeisel_exp)
    storage.mode(zeisel_exp) <- "numeric"
    # Extract the annotation data
    zeisel_annot <- all_dat[seq_len((expDataStarts - 2)), 2:dim(all_dat)[2]]
    zeisel_annot <- rbind(zeisel_annot, tissue = colnames(zeisel_annot))
    rownames(zeisel_annot) <- zeisel_annot[, 1]
    zeisel_annot <- zeisel_annot[, -1]
    z2_annot <- data.frame(
        groupNo = as.numeric(zeisel_annot[1, ]),
        total_mRNA_mol = as.numeric(zeisel_annot[2, ]),
        well = as.numeric(zeisel_annot[3, ]),
        sex = as.numeric(zeisel_annot[4, ]),
        age = as.numeric(zeisel_annot[5, ]),
        diameter = as.numeric(zeisel_annot[6, ]),
        cell_id = as.character(zeisel_annot[7, ]),
        level1class = as.character(zeisel_annot[8, ]),
        level2class = as.character(zeisel_annot[9, ]),
        tissue = as.character(zeisel_annot[10, ]), stringsAsFactors = FALSE
    )
    z2_annot$tissue <- gsub("\\..*", "", as.character(z2_annot$tissue))
    colnames(zeisel_exp) <- z2_annot$cell_id
    rownames(zeisel_exp) <- all_dat[expDataStarts:dim(all_dat)[1], 1]
    rownames(z2_annot) <- z2_annot$cell_id
    return(list(exp = zeisel_exp, annot = z2_annot))
}
