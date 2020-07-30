download.file("goo.gl/r5Y24y",
              destfile="expression_mRNA_17-Aug-2014.txt")

path = "expression_mRNA_17-Aug-2014.txt"

library(EWCE)
cortex_mrna  = load.linnarsson.sct.data(path)

#####
generate.celltype.data <- function(exp,annot,numCores=1){
    library(Matrix)
    library(future)

    # Convert to data.table and inner join so it's clear what expression is from which cell type (this can be very slow with base R so use data.table)
    # - guide to data.table merge functions is here: https://gist.github.com/nacnudus/ef3b22b79164bbf9c0ebafbf558f22a0
    if(class(exp)=="matrix"){
        exp_dt = data.table(exp,keep.rownames = TRUE)
        exp_sparse = Matrix(exp)
    }else{
        stop("Expected exp to be a normal matrix... need to write functions for other instances")
    }

    # Set parameters to enable scTransform to use multiple cores
    future::plan(strategy = 'multicore', workers = numCores)
    options(future.globals.maxSize = 10 * 1024 ^ 3)

    # Variable transform the data using scTransform
    #devtools::install_github(repo = 'ChristophH/sctransform')
    normalized_exp <- sctransform::vst(exp_sparse,return_corrected_umi=TRUE) #$umi_corrected
    normalized_exp_umi = normalized_exp$umi_corrected
    normalized_exp_dt = data.table(as.matrix(normalized_exp_umi),keep.rownames = TRUE) # Need to convert into a normal matrix here.... this will slow the function down  when using large datasets

    annot = data.table(annot)
    exp_long   = data.table::melt(normalized_exp_dt) %>% dplyr::rename(GeneSymbol=rn,cell_id=variable,Exp=value)
    exp_merged = exp_long[annot, on="cell_id", nomatch=0]

    # Use data.table's fast aggregate functionality to get mean expression of each gene, in each celltype
    ct_means = exp_merged[,list(mean=mean(Exp)),by = list(GeneSymbol,level1class)]

    # Calculate total mean expression across all cell types
    ct_totalExp = ct_means[,list(totalExp=sum(mean)),by = list(GeneSymbol)]

    # Merge back together
    ct_meansSummed = ct_means[ct_totalExp, on="GeneSymbol", nomatch=0]

    # Calculate specificity
    ct_specificity = ct_meansSummed[,list(specificity=mean/totalExp),by = list(GeneSymbol,level1class)]

    return(ct_specificity)
}

exp = cortex_mrna$exp
annot = cortex_mrna$annot
ct_specificity = generate.celltype.data(exp,annot,4)
#####



get_ct_means_lvl1 = function(x){
    ct_means = t(aggregate(x,by=list(cell.attrs[["level1class"]]),FUN=mean)[,-1])
    return(t(ct_means))
}

library(Matrix.utils)
ct_means = t(Matrix.utils::aggregate.Matrix(cortex_mrna$exp,groupings=as.factor(cortex_mrna$annot$level1class),FUN=mean))#[,-1])

# Convert it to have a long format... and then aggregate
# https://davetang.org/muse/2016/10/13/using-dplyr-aggregate-r/




exp_lvl5 <- cortex_mrna$exp %>% group_by(Gene,cell_type) %>% summarise(sumUMI=sum(sumUMI)) %>% ungroup()

