#' Create background gene list for multiple species
#'
#' Create background gene list for the
#' intersection/union between multiple species
#'  (\code{gene_list1_species}, \code{gene_list2_species}, and
#'  \code{sctSpecies}), and then filter the gene lists to only include genes
#'  within the background.
#'
#' @inheritParams orthogene::create_background
#' @inheritParams bootstrap_enrichment_test
#'
#' @return Background and gene list.
#'
#' @keywords internal
#' @importFrom orthogene convert_orthologs create_background
create_background_multilist <- function(gene_list1,
                                        gene_list2,
                                        gene_list1_species,
                                        gene_list2_species,
                                        output_species = "human",
                                        bg = NULL,
                                        use_intersect = FALSE,
                                        method = "homologene",
                                        verbose = TRUE) {
    #### If all species are the same, just use all_genes ####
    if (all(c(gene_list1_species, gene_list2_species) == output_species)) {
        gene_map <- orthogene::all_genes(
            species = output_species,
            method = method,
            verbose = verbose
        )
        bg <- gene_map$Gene.Symbol
        return(list(
            bg = bg,
            gene_list1 = gene_list1,
            gene_list2 = gene_list2
        ))
    }
    #### Gene list 1 ####
    #### check gene_list1 args #####
    if (gene_list1_species != output_species) {
        gene_list1 <- orthogene::convert_orthologs(
            gene_df = gene_list1,
            input_species = gene_list1_species,
            output_species = output_species,
            gene_output = "dict"
        )
    }
    bg1 <- orthogene::create_background(
        species1 = gene_list1_species,
        species2 = output_species,
        method = "homologene",
        verbose = verbose
    )

    #### Gene list 2 ####
    #### check gene_list2 args #####
    if (gene_list2_species != output_species) {
        gene_list2 <- orthogene::convert_orthologs(
            gene_df = gene_list2,
            input_species = gene_list2_species,
            output_species = output_species,
            gene_output = "dict"
        )
    }
    if (gene_list1_species == gene_list2_species) {
        bg <- bg2 <- bg1
    } else {
        #### Create bg2 only if different from bg1 ####
        bg2 <- orthogene::create_background(
            species1 = gene_list2_species,
            species2 = output_species,
            method = "homologene",
            verbose = verbose
        )
        #### Create union background ####
        if (use_intersect) {
            bg <- intersect(bg1, bg2)
        } else {
            bg <- union(bg1, bg2)
        }
    }
    #### Report ####
    messager("Using", if (use_intersect) "intersect" else "union",
        "between background gene lists:",
        formatC(length(bg), big.mark = ","), "genes.",
        v = verbose
    )
    return(list(
        bg = bg,
        gene_list1 = gene_list1,
        gene_list2 = gene_list2
    ))
}
