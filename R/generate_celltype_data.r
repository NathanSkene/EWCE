#' Generate CellTypeData (CTD) file
#'
#' \code{generate_celltype_data} takes gene expression data and
#' cell type annotations and creates CellTypeData (CTD) files which
#' contain matrices of mean expression and specificity per cell type.
#'
#' @param exp Numerical matrix with row for each gene and column for each cell.
#' Row names are gene symbols. Column names are cell IDs which can be
#' cross referenced against the annot data frame.
#' @param annotLevels List with arrays of strings containing the cell type
#' names associated with each column in \code{exp}.
#' @param groupName A human readable name for referring to the dataset being
#' @param no_cores Number of cores that should be used to speedup the
#' computation.
#' \emph{NOTE}: Use \code{no_cores=1} when using this package in windows system.
#' @param savePath Directory where the CTD file should be saved.
#' @param file_prefix Prefix to add to saved CTD file name.
#' @param as_sparse Convert \code{exp} to a sparse \code{Matrix}.
#' @param as_DelayedArray Convert \code{exp} to \code{DelayedArray}.
#' @param convert_orths If \code{input_species!=output_species} and
#' \code{convert_orths=TRUE}, will drop genes without
#' 1:1 \code{output_species} orthologs and then convert \code{exp} gene names
#' to those of \code{output_species}.
#' @param input_species The species that the \code{exp} dataset comes from.
#' See \link[EWCE]{list_species} for all available species.
#' @param output_species Species to convert \code{exp} to
#' (Default: "human").
#' See \link[EWCE]{list_species} for all available species.
#' @param force_new_file If a file of the same name as the one
#' being created already exists, overwrite it.
#' @param specificity_quantiles Compute specificity quantiles.
#' Recommended to set to \code{TRUE}.
#' @param dendrograms Add dendrogram plots
#' @param return_ctd Return the CTD object in a list along with the file name,
#' instead of just the file name.
#' @param normSpec Boolean indicating whether specificity data should be
#' transformed to a normal distribution by cell type, giving equivalent scores
#' across all cell types.
#' @param verbose Print messages.
#' @inheritParams bin_specificity_into_quantiles
#' @inheritParams prep_dendro
#' @inheritParams orthogene::convert_orthologs
#' @inheritDotParams orthogene::convert_orthologs
#'
#' @return File names for the saved CellTypeData (CTD) files.
#' @examples
#' # Load the single cell data
#' cortex_mrna <- ewceData::cortex_mrna()
#' # Use only a subset to keep the example quick
#' expData <- cortex_mrna$exp[1:100, ]
#' l1 <- cortex_mrna$annot$level1class
#' l2 <- cortex_mrna$annot$level2class
#' annotLevels <- list(l1 = l1, l2 = l2)
#' fNames_ALLCELLS <- EWCE::generate_celltype_data(
#'     exp = expData,
#'     annotLevels = annotLevels,
#'     groupName = "allKImouse"
#' )
#' @export
#' @importFrom parallel mclapply
#' @importFrom Matrix Matrix
#' @importFrom RNOmni RankNorm
#' @importFrom methods as is
#' @importFrom orthogene convert_orthologs
generate_celltype_data <- function(exp,
                                   annotLevels,
                                   groupName,
                                   no_cores = 1,
                                   savePath = tempdir(),
                                   file_prefix = "ctd",
                                   as_sparse = TRUE,
                                   as_DelayedArray = FALSE,
                                   normSpec = FALSE,
                                   convert_orths = FALSE,
                                   input_species = "mouse",
                                   output_species = "human",
                                   non121_strategy = "drop_both_species",
                                   method = "homologene",
                                   force_new_file = TRUE,
                                   specificity_quantiles = TRUE,
                                   numberOfBins = 40,
                                   dendrograms = TRUE,
                                   return_ctd = FALSE,
                                   verbose = TRUE,
                                   ...) {
    #### Check if input is an SCE or SE object ####
    res_sce <- check_sce(exp)
    exp <- res_sce$exp
    SE_obj <- res_sce$SE_obj
    #### Check if file already exists ####
    fNames <- sprintf("%s/%s_%s.rda", savePath, file_prefix, groupName)
    if (file.exists(fNames) & (!force_new_file)) {
        messager("+ Pre-existing file of the same name detected.",
            "Use `force_new_file=TRUE` to overwrite this file.",
            v = verbose
        )
        messager("+ Returning pre-existing file path.", v = verbose)
        return(fNames)
    }
    #### Assign cores #####
    core_allocation <- assign_cores(
        worker_cores = no_cores,
        verbose = verbose
    )
    no_cores <- core_allocation$worker_cores
    #### Check NAs ####
    check_nas(exp = exp)
    #### Check annotLevels ####
    check_annotLevels(
        annotLevels = annotLevels,
        exp = exp
    )
    #### Check groupName ####
    check_group_name(groupName = groupName)
    #### convert orthologs ####
    if ((input_species != output_species) && convert_orths) {
        exp <- orthogene::convert_orthologs(
            gene_df = exp,
            input_species = input_species,
            output_species = output_species,
            method = method,
            non121_strategy = non121_strategy,
            as_sparse = as_sparse,
            as_DelayedArray = as_DelayedArray,
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
    ### Avoid conflicts with parallelized block-wise operations ####
    if (is_delayed_array(exp)) {
        no_cores <- 1
    }
    #### Initialize ctd ####
    ctd <- list()
    for (i in seq_len(length(annotLevels))) {
        ctd[[length(ctd) + 1]] <- list(annot = annotLevels[[i]])
    }
    #### Calculate mean_exp and specificity ####
    messager("+ Calculating normalized mean expression.", v = verbose)
    ctd2 <- parallel::mclapply(ctd, calculate_meanexp_for_level,
        exp,
        mc.cores = no_cores
    )
    messager("+ Calculating normalized specificity.", v = verbose)
    ctd3 <- parallel::mclapply(ctd2, calculate_specificity_for_level,
        mc.cores = no_cores
    )
    ctd <- ctd3
    remove(ctd2, ctd3)
    #### Apply rank norm transformation if hyp param set ####
    if (isTRUE(normSpec)) {
        ctd <- lapply(ctd, rNorm)
    }
    #### Calculate specificity quantiles ####
    if (specificity_quantiles) {
        ctd <- lapply(ctd,
            bin_specificity_into_quantiles,
            numberOfBins = numberOfBins
        )
    }
    #### Add dendrograms ####
    if (dendrograms) {
        ctd <- lapply(ctd, prep_dendro)
    }
    #### Save results ####
    messager("+ Saving results ==> ", fNames, v = verbose)
    dir.create(dirname(fNames), showWarnings = FALSE, recursive = TRUE)
    save(ctd, file = fNames)
    #### Return ####
    if (return_ctd) {
        messager("+ Returning list of CTD file name, and the CTD itself.",
            v = verbose
        )
        return(list(
            fNames = fNames,
            ctd = ctd
        ))
    } else {
        return(fNames)
    }
}
