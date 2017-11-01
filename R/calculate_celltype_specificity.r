#' Calculate celltype specificity
#'
#' \code{calculate_celltype_specificity} takes expression & annotation data and determines
#' for each gene, the proportion of expression found in each celltype
#'
#' @param exp A matrix containing the expression data, wherein rows are genes and columns are cells
#' @param annot A list of annotation vectors, where each vector has a name for each cell
#' @param trim Value determining how the trimmed mean is calculated (range: 0 to 0.5)
#' @param thresh Expression threshold value below which a gene is dropped
#' @return A list array, in which each cell contains two data frames:
#' \itemize{
#'   \item \code{mean_exp}: stores the average value for each subcell type found across cells
#'   \item \code{specificity}: proportion of expression for each gene in each cell type
#' }
#' @export
#' @examples
#' # Load the single cell data
#' data("cortex_mrna")
#' ctd = calculate_celltype_specificity(cortex_mrna$exp,cortex_mrna$annot$level2class)
calculate_celltype_specificity <- function(exp,annot,thresh=0,trim=0){
    # Check that annot vectors have length equal to the number of cells
    for(i in 1:length(annot)){
        if(!length(annot)==(dim(exp)[2])){stop(sprintf("Length of annotation vector %s is not equal to the number of cells",i))}
    }

    # Convert characters to numbers
    expr2 = matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),ncol=ncol(exp))
    rownames(expr2) = rownames(exp)
    colnames(expr2) = colnames(exp)

    celltype_data = list()
    #for(i in 1:length(annot)){celltype_data[[i]]=list()}

    # Get mean values for subcell types
    count = 0
    for(sct in unique(annot)){
        count = count+1
        sub_expr = expr2[,annot==sct]
        sct_expr = data.frame(temp=apply(sub_expr,1,mean,trim=trim))
        colnames(sct_expr)=sct
        if(count==1){
            all_scts = sct_expr
        }else{
            all_scts = cbind(all_scts,sct_expr)
        }
    }

    # Drop SCT genes with very low levels of expression
    keepGenes = rownames(all_scts)[apply(all_scts,1,max)>thresh]
    all_scts = all_scts[keepGenes,]

    # Get cell type proportion for cell type level annotations
    cTs = unique(annot)
    geneList= unique(rownames(all_scts))
    count = 0
    for(gs in geneList){
        count=count+1

        exp1 = unlist(all_scts[gs,])

        exp2 = exp1/sum(exp1) #<- This is just used to create the vector
        for(ct in unique(annot)){
            #exp2[ct] = exp1[ct]/sum(exp1[a0[[ct]]$annot])
            if(is.nan(exp2[ct])){exp2[ct]=0}
            if(is.infinite(exp2[ct])){exp2[ct]=0}
        }

        exp5 = t(data.frame(exp2))
        rownames(exp5) = gs
        if(count==1){
            cell_dists = exp5
        }else{
            cell_dists = rbind(cell_dists,exp5)
        }
    }

    celltype_data$mean_exp = all_scts
    celltype_data$specificity = cell_dists

    return(celltype_data)
}
