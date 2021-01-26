#' merge_two_expfiles
#'
#' \code{merge_two_expfiles} Used to combine two single cell type datasets
#'
#' @param exp1 Numerical expression matrix for dataset1 with row for each gene and column for each cell. Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be cross referenced against the annot data frame.
#' @param exp2 Numerical expression matrix for dataset2 with row for each gene and column for each cell. Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be cross referenced against the annot data frame.
#' @param annot1 Annotation data frame for dataset1 which contains three columns at least: cell_id, level1class and level2class
#' @param annot2 Annotation data frame for dataset2 which contains three columns at least: cell_id, level1class and level2class
#' @param name1 Name used to refer to dataset 1. Leave blank if it's already a merged dataset.
#' @param name2 Name used to refer to dataset 2. Leave blank if it's already a merged dataset.
#' @return List containing merged exp and annot
#' @examples
#' # See the EWCE vignette
#' @export
merge_two_expfiles <- function(exp1,exp2,annot1,annot2,name1="",name2=""){
    #genes_intersect   = intersect(rownames(exp1),rownames(exp2))
    #genes_exp1_unique = setdiff(rownames(exp1),rownames(exp2))
    #genes_exp2_unique = setdiff(rownames(exp2),rownames(exp1))
    all_genes = unique(c(rownames(exp2),rownames(exp1)))

    # Check for correct column headers in annot data tables
    if(sum(c("cellid") %in% colnames(annot1))==1){colnames(annot1)[colnames(annot1)=="cellid"]="cell_id"}
    if(sum(c("cellid") %in% colnames(annot2))==1){colnames(annot2)[colnames(annot2)=="cellid"]="cell_id"}
    if(!sum(c("cell_id","level1class","level2class") %in% colnames(annot1))==3){stop("ERROR: annot1 doesn't have either cell_id, level1class or level2class columns")}
    if(!sum(c("cell_id","level1class","level2class") %in% colnames(annot2))==3){stop("ERROR: annot2 doesn't have either cell_id, level1class or level2class columns")}
    
    # If one of the exp matrices is really a matrix (not a dataframe) the code won't work... so force conversion
    #if(class(exp1)=="matrix"){exp1=data.frame(exp1);colnames(exp1)=gsub("X","",colnames(exp1))}
    #if(class(exp2)=="matrix"){exp2=data.frame(exp2);colnames(exp2)=gsub("X","",colnames(exp2))}
    if(class(exp1)=="matrix"){nn=colnames(exp1);rr=rownames(exp1);exp1=data.frame(exp1,stringsAsFactors = FALSE);colnames(exp1)=nn;rownames(exp1)=rr;}
    if(class(exp2)=="matrix"){nn=colnames(exp2);rr=rownames(exp2);exp2=data.frame(exp2,stringsAsFactors = FALSE);colnames(exp2)=nn;rownames(exp2)=rr;}


    # Merge the expression matrices, setting undetected genes to 0
    exp1b = exp1[all_genes,]
    exp1b[is.na(exp1b)]=0
    rownames(exp1b)=all_genes
    exp2b = exp2[all_genes,]
    exp2b[is.na(exp2b)]=0
    rownames(exp2b)=all_genes
    exp =cbind(exp1b,exp2b)

    # Ensure important annotation columns are not stored as factors
    annot1$cell_id     = as.character(annot1$cell_id)
    annot2$cell_id     = as.character(annot2$cell_id)
    annot1$level1class = as.character(annot1$level1class)
    annot1$level2class = as.character(annot1$level2class)
    if(is.null(annot1$dataset_name)){
        annot1$dataset_name = name1
    }
    annot2$level1class = as.character(annot2$level1class)
    annot2$level2class = as.character(annot2$level2class)
    if(is.null(annot2$dataset_name)){
        annot2$dataset_name = name2
    }
    keepTISSUE=FALSE
    #if(("tissue" %in% colnames(annot1)) & ("tissue" %in% colnames(annot2))){
    if(all(("tissue" %in% colnames(annot1)) & ("tissue" %in% colnames(annot2)))){
        keepTISSUE=TRUE    
        annot1$tissue     = as.character(annot1$tissue)
        annot2$tissue     = as.character(annot2$tissue)
    }

    # Setup new annotation dataframe
    cell_id = c(annot1$cell_id,annot2$cell_id)
    level1class = c(annot1$level1class,annot2$level1class)
    level2class = c(annot1$level2class,annot2$level2class)
    dataset_name = c(as.character(annot1$dataset_name),as.character(annot2$dataset_name))
    if(!keepTISSUE){
        annot = data.frame(cell_id=cell_id,level1class=level1class,level2class=level2class,dataset_name=dataset_name)    
    }else{
        tissue = c(annot1$tissue,annot2$tissue)
        annot = data.frame(cell_id=cell_id,level1class=level1class,level2class=level2class,dataset_name=dataset_name,tissue=tissue)
    }
    

    # Drop expression data without annotation data
    numMissingAnnot = dim(exp)[2] - sum(annot$cell_id %in% colnames(exp))
    if(numMissingAnnot>0){
        print(sprintf("Warning: %s cells are missing annotation data and have been dropped",numMissingAnnot))
        exp = exp[,as.character(annot$cell_id)]
    }

    return(list(exp=exp,annot=annot))
}
