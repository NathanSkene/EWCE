#' Drop uninformative genes
#'
#'
#' \code{drop_uninformative_genes} drops uninformative genes in order to reduce
#' compute time and noise in subsequent steps. It achieves this through several
#' steps, each of which are optional:
#' \itemize{
#' \item{Drop non-1:1 orthologs:\cr}{Removes genes that don't have 1:1 orthologs
#' with the \code{output_species} ("human" by default).}
#' \item{Drop non-varying genes:\cr}{Removes genes that don't vary across cells
#' based on variance deciles.}
#' \item{Drop non-differentially expressed genes (DEGs):\cr}{
#' Removes genes that are not significantly differentially
#' expressed across cell-types (multiple DEG methods available).}
#' }
#'
#' @param exp Expression matrix with gene names as rownames.
#' @param level2annot Array of cell types, with each sequentially corresponding
#' a column in the expression matrix.
#' @param mtc_method Multiple-testing correction method used by DGE step.
#' See \link[stats]{p.adjust} for more details.
#' @param return_sce Whether to return the filtered results
#' as an expression matrix or a \pkg{SingleCellExperiment}.

#' @param adj_pval_thresh Minimum differential expression significance
#' that a gene must demonstrate across \code{level2annot} (i.e. cell types).
#' @param input_species Which species the gene names in \code{exp} come from.
#' See \link[EWCE]{list_species} for all available species.
#' @param output_species Which species' genes names to convert \code{exp} to.
#' See \link[EWCE]{list_species} for all available species.
#' @param as_sparse Convert \code{exp} to sparse matrix.
#' @param as_DelayedArray Convert \code{exp} to \code{DelayedArray}
#'  for scalable processing.
#' @param no_cores Number of cores to parallelise across.
#' Set to \code{NULL} to automatically optimise.
#' @param verbose Print messages.
#' @inheritParams orthogene::convert_orthologs
#' @inheritParams generate_celltype_data
#' @inheritDotParams orthogene::convert_orthologs
#'
#' @return exp Expression matrix with gene names as row names.
#' @examples
#' cortex_mrna <- ewceData::cortex_mrna()
#' # Use only a subset of genes to keep the example quick
#' cortex_mrna$exp <- cortex_mrna$exp[1:300, ]
#'
#' ## Convert orthologs at the same time
#' exp2_orth <- drop_uninformative_genes(
#'     exp = cortex_mrna$exp,
#'     level2annot = cortex_mrna$annot$level2class,
#'     input_species = "mouse"
#' )
#' @export
#' @importFrom Matrix rowSums colSums
#' @importFrom methods is as
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom DelayedArray DelayedArray
#' @importFrom stats p.adjust
#' @importFrom orthogene convert_orthologs
drop_uninformative_genes <- function(exp,
                                     level2annot,
                                     #dge_method = "limma",
                                     #dge_test = "LRT",
                                     mtc_method = "BH",
                                     #min_variance_decile = NULL,
                                     adj_pval_thresh = 0.00001,
                                     convert_orths = FALSE,
                                     input_species = NULL,
                                     output_species = "human",
                                     non121_strategy = "drop_both_species",
                                     method = "homologene",
                                     as_sparse = TRUE,
                                     as_DelayedArray = FALSE,
                                     return_sce = FALSE,
                                     no_cores = 1,
                                     verbose = TRUE,
                                     ...) {
    
    ##### Extra arguments to be implemented after benchmarking is done #### 
  #' @param dge_method Which method to use for the Differential Gene Expression
  #' (DGE) step.\cr
  #' Options:
  #' \itemize{
  #' \item{"limma": }{Uses
  #' \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{
  #' limma}.}
  #' \item{"deseq2": }{Uses
  #'  \href{https://bioconductor.org/packages/release/bioc/html/DESeq2.html}{
  #' DESeq2}.}
  #' \item{"mast": }{Uses
  #' \href{https://www.bioconductor.org/packages/release/bioc/html/MAST.html}{
  #' MAST}.}
  #' }
  #' @param dge_test \code{test} argument passed to DGE function.
  #' Only used when \code{dge_method="deqseq2"}.
  #' @param min_variance_decile If \code{min_variance_decile!=NULL},
  #'  calculates the variance of the mean gene expression
  #'   across `level2annot` (i.e. cell types),
  #' and then removes any genes that are below \code{min_variance_decile}
  #'  (on a 0-1 scale).
    dge_method = "limma"
    dge_test = "LRT"
    min_variance_decile = NULL
    #### End extra arguments ####
    
    ### Avoid confusing Biocheck
    padj <- NULL
    #### Check if input is an SCE or SE object ####
    res_sce <- check_sce(exp)
    exp <- res_sce$exp
    SE_obj <- res_sce$SE_obj
    metadata <- res_sce$metadata
    messager("Check", dim(exp), v = verbose)
    #### Check species ####
    input_species <- check_species(
        genelistSpecies = "NULL",
        sctSpecies = input_species,
        verbose = verbose
    )$sctSpecies
    #### Check DGE method ####
    dge_method <- if (is.null(dge_method)) "" else dge_method
    #### Assign cores #####
    core_allocation <- assign_cores(
        worker_cores = no_cores,
        verbose = verbose
    )
    #### convert orthologs ####
    if ((input_species != output_species) && convert_orths) {
        exp <- orthogene::convert_orthologs(
            gene_df = exp,
            input_species = input_species,
            output_species = output_species,
            non121_strategy = non121_strategy,
            method = method,
            verbose = verbose,
            ...
        )
    }
    #### Convert to sparse matrix ####
    exp <- to_sparse_matrix(
        exp = exp,
        as_sparse = as_sparse,
        verbose = verbose
    )
    #### Convert to DelayedArray ####
    exp <- to_delayed_array(
        exp = exp,
        as_DelayedArray = as_DelayedArray,
        verbose = verbose
    )

    ### Remove non-expressed genes ####
    exp <- drop_nonexpressed_genes(
        exp = exp,
        verbose = verbose
    )
    ### Remove low-quality cells ####
    exp_annotLevels <- drop_nonexpressed_cells(
        exp = exp,
        annotLevels = list(lvl2 = level2annot),
        verbose = verbose
    )
    exp <- exp_annotLevels$exp
    level2annot <- exp_annotLevels$annotLevels[[1]]

    #### Simple variance ####
    if (!is.null(min_variance_decile)) {
        # Use variance of mean gene expression across cell types
        # as a fast and simple way to select genes
        exp <- filter_variance_quantiles(
            exp = exp,
            n_quantiles = 10,
            min_variance_quantile = min_variance_decile,
            verbose = verbose
        )
    }
    #### Run DGE ####
    start <- Sys.time()

    #### limma ####
    if (any(tolower(dge_method) == "limma")) {
        limma_res <- run_limma(
            exp = exp,
            level2annot = level2annot,
            mtc_method = mtc_method,
            verbose = verbose,
            ...
        )
        keep_genes <- rownames(exp)[limma_res$q < adj_pval_thresh]
        report_dge(
            exp = exp,
            keep_genes = keep_genes,
            adj_pval_thresh = adj_pval_thresh,
            verbose = verbose
        )
        # Filter original exp
        exp <- exp[keep_genes, ]
    }

    #### DESeq2 ####
    if (tolower(dge_method) == "deseq2") {
        deseq2_res <- run_deseq2(
            exp = exp,
            level2annot = level2annot,
            test = dge_test,
            no_cores = no_cores,
            verbose = verbose,
            ...
        )
        keep_genes <- rownames(subset(deseq2_res, padj < adj_pval_thresh))
        report_dge(
            exp = exp,
            keep_genes = keep_genes,
            adj_pval_thresh = adj_pval_thresh,
            verbose = verbose
        )
        # Filter original exp
        exp <- exp[keep_genes, ]
    }

    #### MAST ####
    if (tolower(dge_method) == "mast") {
        mast_res <- run_mast(
            exp = exp,
            level2annot,
            test = dge_test,
            mtc_method = mtc_method,
            no_cores = no_cores,
            ...
        )
        keep_genes <- unique(subset(mast_res, q < adj_pval_thresh)$primerid)
        report_dge(
            exp = exp,
            keep_genes = keep_genes,
            adj_pval_thresh = adj_pval_thresh,
            verbose = verbose
        )
        # Filter original exp
        exp <- exp[keep_genes, ]
    }
    end <- Sys.time()
    print(end - start) ## End DGE

    #### Return results ####
    if (return_sce) {
        new_sce <- SingleCellExperiment::SingleCellExperiment(
            list(counts = exp),
            colData = metadata[colnames(exp), ]
        )
        return(new_sce)
    } else {
        return(exp)
    }
}



### OG drop.uninformative.genes ####
drop.uninformative.genes <- function(exp,
                                     level2annot,
                                     ...) {
    .Deprecated("filter_nonorthologs")
    exp <- drop_uninformative_genes(
        exp = exp,
        level2annot = level2annot,
        ...
    )
    # level2annot = as.character(level2annot)
    # summed = apply(exp,1,sum)
    # exp = exp[summed!=0,]
    # mod  = model.matrix(~level2annot)
    # fit = limma::lmFit(exp,mod)
    # eb = limma::eBayes(fit)
    # pF = stats::p.adjust(eb$F.p.value,method="BH")
    # exp = exp[pF<0.00001,]
    # return(exp)
}
