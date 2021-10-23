#' Drop uninformative genes
#'
#' \code{drop_uninformative_genes} first removes genes
#' that do not have 1:1 orthologs with humans.
#'
#' \code{drop_uninformative_genes} then drops genes from an
#' SCT expression matrix
#' if they do not significantly vary between any cell types.
#' Makes this decision based on use of an ANOVA (implemented with \code{limma}).
#' If the F-statistic for variation amongst type2 annotations
#' is less than a strict p-threshold, then the gene is dropped.
#'
#' A very fast alternative to DGE methods is filtering by
#' \code{min_variance_decile},
#' which selects only genes with the top variance deciles.
#'
#' @param exp Expression matrix with gene names as rownames.
#' @param level2annot Array of cell types, with each sequentially corresponding
#' a column in the expression matrix
#' @param DGE_method Which method to use for the Differential Gene Expression
#' (DGE) step.
#' @param return_sce Whether to return the filtered results
#' as an expression matrix or a \pkg{SingleCellExperiment}.
#' @param min_variance_decile If \code{min_variance_decile!=NULL},
#'  calculates the variance of the mean gene expression
#'   across `level2annot` (i.e. cell types),
#' and then removes any genes that are below \code{min_variance_decile}
#'  (on a 0-1 scale).
#' @param adj_pval_thresh Minimum differential expression significance
#' that a gene must demonstrate across \code{level2annot} (i.e. cell types).
#' @param input_species Which species the gene names in \code{exp} come from.
#' @param output_species Which species' genes names to convert \code{exp} to.
#' @param as_sparse Convert \code{exp} to sparse matrix.
#' @param as_DelayedArray Convert \code{exp} to \code{DelayedArray}
#'  for scalable processing.
#' @param no_cores Number of cores to parallelise across.
#' Set to \code{NULL} to automatically optimise.
#' @param verbose Print messages.
#' @param ... Additional arguments to be passed to the selected DGE method.
#' @inheritParams orthogene::convert_orthologs
#' @inheritParams generate_celltype_data
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
                                     DGE_method = "limma",
                                     min_variance_decile = NULL,
                                     adj_pval_thresh = 0.00001,
                                     convert_orths = FALSE,
                                     input_species = NULL,
                                     output_species = "human",
                                     non121_strategy = "drop_both_species",
                                     as_sparse = TRUE,
                                     as_DelayedArray = FALSE,
                                     return_sce = FALSE,
                                     no_cores = 1,
                                     verbose = TRUE,
                                     ...) {
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
    DGE_method <- if (is.null(DGE_method)) "" else DGE_method
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
            verbose = verbose,
            ...
        )
    }
    #### Convert to sparse matrix ####
    exp <- to_sparse_matrix(
        exp = exp,
        as_sparse = as_sparse
    )
    #### Convert to DelayedArray ####
    exp <- to_delayed_array(
        exp = exp,
        as_DelayedArray = as_DelayedArray
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
            level2annot = level2annot,
            n_quantiles = 10,
            min_variance_quantile = min_variance_decile,
            verbose = verbose
        )
    }
    #### Run DGE ####
    start <- Sys.time()
    #### Limma ####
    # Modified original method
    if (any(tolower(DGE_method) == "limma")) {
        eb <- run_limma(
            exp = exp,
            level2annot = level2annot,
            verbose = verbose,
            ...
        )
        pF <- stats::p.adjust(eb$F.p.value, method = "BH")
        keep_genes <- pF < adj_pval_thresh
        messager(paste(
            nrow(exp) - sum(keep_genes), "/", nrow(exp),
            "genes dropped @ DGE adj_pval_thresh <", adj_pval_thresh
        ), v = verbose)
        # Filter original exp
        exp <- exp[keep_genes, ]
    }

    #### DESeq2 ####
    if (tolower(DGE_method) == "deseq2") {
        dds_res <- run_deseq2(
            exp = exp,
            level2annot = level2annot,
            verbose = verbose,
            ...
        )
        dds_res <- subset(dds_res, padj < adj_pval_thresh)
        messager(paste(
            nrow(exp) - nrow(dds_res), "/", nrow(exp),
            "genes dropped @ DGE adj_pval_thresh <", adj_pval_thresh
        ), v = verbose)
        # Filter original exp
        exp <- exp[row.names(dds_res), ]
    }

    #### glmGamPoi ####
    ## Removing this option for now until we can figure out
    ## how to pass the Travis CI checks,
    ## which are failing when installing the deps for glmGamPo (hdf5).
    # if(tolower(DGE_method)=="glmgampoi"){
    #     # Best for large datasets that can't fit into memory.
    #     messager("DGE:: glmGamPoi...",v=verbose)
    #     sce_de <- run_glmGamPoi_DE(sce,
    #                                level2annot=level2annot,
    #                                adj_pval_thresh=adj_pval_thresh,
    #                                return_as_SCE=T,
    #                                # sce_save_dir=sce_save_dir,
    #                                verbose=verbose)
    #                                # ...)
    #     messager(paste(
    #         nrow(sce)-nrow(sce_de),"/",nrow(sce),
    #         "genes dropped @ DGE adj_pval_thresh <",adj_pval_thresh),
    #         v=verbose)
    #     sce <- sce_de
    # }
    # Report time elapsed
    end <- Sys.time()
    print(end - start)
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
