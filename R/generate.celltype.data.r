#' generate.celltype.data
#'
#' \code{generate.celltype.data} Takes expression & cell type annotations and creates celltype_data files which contain the mean and specificity matrices
#'
#' @param exp Numerical matrix with row for each gene and column for each cell. Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be cross referenced against the annot data frame.
#' @param level1class Array of strings containing the level 1 cell type names associated with each column in exp
#' @param level2class Array of strings containing the level 2 cell type names associated with each column in exp
#' @param groupName A human readable name for refering to the dataset being loaded
#' @param trim Value determining how the trimmed mean is calculated (range: 0 to 0.5, default=0). We do not recommend changing from the default.
#' @param thresh Expression threshold value below which a gene is dropped (default=0). We do not recommend changing from the default.
#' @return Filenames for the saved celltype_data files
#' @examples
#' # Load the single cell data
#' data("cortex_mrna")
#' expData = cortex_mrna$exp
#' l1=cortex_mrna$annot$level1class
#' l2=cortex_mrna$annot$level2class
#' fNames_ALLCELLS = generate.celltype.data(exp=expData,l1,l2,"allKImouse",thresh=0,trim=0)
#' @export
generate.celltype.data <- function(exp,level1class,level2class,groupName,thresh=0,trim=0){
    # First, check the number of annotations equals the number of columns in the expression data
    if(length(level1class)!=length(level2class)){stop("Error: length of level1class does not equal that of level2class")}
    if(length(level1class)!=dim(exp)[2]){stop("Error: length of level1class must equal the number of columns in exp matrix")}

    # Check group name
    if(is.null(groupName)){stop("ERROR: groupName must be set. groupName is used to label the files created by this function.")}
    if(groupName==""){stop("ERROR: groupName must be set. groupName is used to label the files created by this function.")}
    
    fNames=rep("",2)
    ctd = list()
    for(lev in 1:2){
        annot_levels = list()
        if(lev==1){   annot_levels = level1class   }
        if(lev==2){   annot_levels = level2class   }
        ctd[[lev]] = calculate_celltype_specificity(exp=exp,annot=annot_levels,thresh=thresh,trim=trim)
        ctd[[lev]]$annot = annot_levels
        acs = ctd[[lev]]$mean_exp
        ctd[[lev]]$mean_exp = acs[,order(colnames(acs))]
        ctd[[lev]]$specificity = ctd[[lev]]$specificity[,order(colnames(acs))]
    }
    if(thresh!=0 | trim!=0){
        fNames=sprintf("CellTypeData_%s_thresh%s_trim%s.rda",groupName)
    }else{
        fNames=sprintf("CellTypeData_%s.rda",groupName)
    }
    save(ctd,file=fNames)
    return(fNames)
}
