#' generate.celltype.data
#'
#' \code{generate.celltype.data} Takes expression & cell type annotations and creates celltype_data files which contain the mean and specificity matrices
#'
#' @param exp Numerical matrix with row for each gene and column for each cell. Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be cross referenced against the annot data frame.
#' @param annotLevels List with arrays of strings containing the cell type names associated with each column in exp
#' @param groupName A human readable name for refering to the dataset being loaded
#' @return Filenames for the saved celltype_data files
#' @examples
#' # Load the single cell data
#' data("cortex_mrna")
#' expData = cortex_mrna$exp
#' l1=cortex_mrna$annot$level1class
#' l2=cortex_mrna$annot$level2class
#' annotLevels = list(l1=l1,l2=l2)
#' fNames_ALLCELLS = generate.celltype.data(exp=expData,annotLevels,"allKImouse")
#' @export
#' @import parallel
generate.celltype.data <- function(exp,annotLevels,groupName){
    require("parallel")
    # Calculate the number of cores
    no_cores <- detectCores()
    cl <- makeCluster(no_cores)

    # First, check the number of annotations equals the number of columns in the expression data
    lapply(annotLevels,test <- function(x,exp){if(length(x)!=dim(exp)[2]){stop("Error: length of all annotation levels must equal the number of columns in exp matrix")}},exp)

    # Check group name
    if(is.null(groupName)){stop("ERROR: groupName must be set. groupName is used to label the files created by this function.")}
    if(groupName==""){stop("ERROR: groupName must be set. groupName is used to label the files created by this function.")}

    # Convert characters to numbers
    exp2<-suppressWarnings(apply(exp,2,function(x) {storage.mode(x) <- 'double'; x}))

    ctd = list()
    for(i in 1:length(annotLevels)){ctd[[length(ctd)+1]] = list(annot=annotLevels[[i]])}

    aggregate.over.celltypes <- function(rowOfMeans,celltypes){
        exp_out = as.matrix(data.frame(aggregate(rowOfMeans,by=list(celltypes),FUN=mean)))
        rownames(exp_out) = exp_out[,"Group.1"]
        exp_out = exp_out[,2]
        exp_out2 = as.numeric(exp_out)
        names(exp_out2) = names(exp_out)
        return(exp_out2)
    }
    calculate.meanexp.for.level <- function(ctd_oneLevel,expMatrix){
        if(dim(expMatrix)[2]==length(unique(ctd_oneLevel$annot))){
            print(dim(expMatrix)[2])
            print(length(ctd_oneLevel$annot))
            if(sum(!colnames(expMatrix)==ctd_oneLevel$annot)!=0){
                stop("There are an equal number of celltypes in expMatrix and ctd_oneLevel but the names do not match")
            }
            ctd_oneLevel$mean_exp = expMatrix
        }else{
            mean_exp = apply(expMatrix,1,aggregate.over.celltypes,ctd_oneLevel$annot)
            ctd_oneLevel$mean_exp = t(mean_exp)
        }
        return(ctd_oneLevel)
    }
    calculate.specificity.for.level <- function(ctd_oneLevel){
        ctd_oneLevel$specificity = ctd_oneLevel$mean_exp/(apply(ctd_oneLevel$mean_exp,1,sum)+0.000000000001)
        return(ctd_oneLevel)
    }
    ctd2 = mclapply(ctd,calculate.meanexp.for.level,exp2)
    ctd3 = mclapply(ctd2,calculate.specificity.for.level)
    ctd=ctd3
    stopCluster(cl)

    fNames=sprintf("CellTypeData_%s.rda",groupName)
    save(ctd,file=fNames)
    return(fNames)
}
