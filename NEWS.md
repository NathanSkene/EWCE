# EWCE 1.15.1

## Bug fixes

* Remove print statements, clean up CTD creation steps.

# EWCE 1.13.1

## Bug fixes

* Making `ewce_plot` functionality more clear - will no fail if 
`make_dendro=TRUE` and `ggdendro` is not installed or `CTD` is not provided 
rather than issuing a warning.

# EWCE 1.11.5

## Bug fixes

* `check_bootstrap_args`:
  - annotation level for checking cell type names was previous hardcoded to 
  level 1, now updated to match user input for annotation level.

# EWCE 1.11.4

## Bug fixes

* `drop_uninformative_genes`:
  - Will now catch cases where expression matrix is a dataframe and convert to a 
  matrix
  - This was causing weird errors, see issue 92: 
  https://github.com/NathanSkene/EWCE/issues/92
  - Fix made in `check_sce()` function.

# EWCE 1.11.2

## New features

* `bootstrap_enrichment_test`:
  - New args: `standardise_sct_data=`, `standardise_hits=`:
  let users have more control over data standardisation steps.
  - `check_ewce_genelist_inputs`: updated accordingly.
  - New arg: `store_gene_data` to avoid hitting memory limits.
  - Modify *test-bootstrap_enrichment_test_2.R* to use test new args.

# EWCE 1.11.1

## Bug fixes

* `ewce_plot()` - Dendrogram not reordering cell types in plot
    - see [issue](https://github.com/NathanSkene/EWCE/issues/88)
    - Occurs when ctd does not contain the plotting info
    - Fixed now and unit test added. 
    - Note that cell type order on the x-axis is based on hierarchical 
      clustering for both plots if `make_dendro = TRUE` for `ewce_plot()`.  

# EWCE 1.9.3

## Bug fixes

* `generate_bootstrap_plots`
    - Use stored `gene_data` object whenever possible.
    - Only show filtered celltypes. 

# EWCE 1.9.2

## Bug fixes

* `generate_bootstrap_plots`
    - Missing available parameters for `check_ewce_genelist_inputs()` call.

# EWCE 1.9.1

## Bug fixes

* `drop_uninformative_genes`
    - Hash out DGE options. Somehow these got re-exposed to users in Bioc>=3.16.

# EWCE 1.7.4

## New features

* `check_species`: 
    - Added new arg `sctSpecies_origin_default`
* `standardise_ctd`:
    - Can now specify `sctSpecies_origin`, which will be added to the metadata.

## Bug fixes

* `standardise_ctd`:
    - When "annot" slot is not provided in the original CTD,
    a new one will be created from matrix columns names instead of assigning `NULL`.

# EWCE 1.7.3

## New features

* `check_ewce_genelist_inputs`/`bootstrap_enrichment_test`
    - New arg: `sctSpecies_origin` lets users clarify that their data originally
    came from mouse even when it is currently formatted as human orthologs. 
    This is necessary for creating the appropriate background gene lists.
* Remove `grDevices` as dep entirely.
* `fix_celltype_names`
    - Add new arg `make_unique` to make this function easily usable for vectors
        where the same celltype appears multiple times.
* `bootstrap_enrichment_test`
    - Return gene-level scores based on adaptation of code from `generate_bootstrap_plots`.
    now stored as a list element named `gene_data` in `data.table` format.
* `generate_bootstrap_plots`
    - Revamp wrap code into reusable subfunctions.
    - Avoid resampling random genes when the gene data is stored in the bootstrap results as `gene_data`.
    It will also tell you which of these options it's using.
    - Save with `ggsave` instead if `grDevices`.
    - Facet by celltype instead of generating tons of separate plots.
    - Let users decide cutoff threshold with new arg `adj_pval_thresh`
    - Now returns a named list with the plots themselves ("plot") and the paths to where they're saved ("paths") rather than just a higher-level directory path in which users had to search for the right files (and didn't ever have access to the ggplot2 objects).
    - Show significance with barplot fill/color instead of asterices. Much easier to see now.
    - Change `savePath` arg to the more accurate `save_dir`. Expose appending *BootstrapPlots* to the user within the argument.
* `generate_bootstrap_plots_for_transcriptome`
    - Change `savePath` arg to the more accurate `save_dir`. Expose appending *BootstrapPlots* to the user within the argument.
    - Save with `ggsave` instead if `grDevices`.
* Standardise `hits` + `hitGenes` arg all to `hits`.
* Update hex:
    - Off load large source image from DALLE to Releases instead of including it within the package.
        
## Bug fixes

* `drop_uninformative_genes` / `generate_celltype_data`
    - Pass `verbose` arg to matrix formatting functions.
* `generate_controlled_bootstrap_geneset`
    - Removed `combinedGenes` arg as it was not being used anywhere within.
* `check_args_for_bootstrap_plot_generation`
    - Removed unused args: `ttSpecies`, `sctSpecies`
* *test-bootstrap_enrichment_test_2.R* 
    - "monkey_ctd" tests seems to be running more smoothly than before (not just getting NAs). 
        This might have to with `orthogene` databases improving.
    - Reassuringly, "godzilla" tests still fail as expected :)
* Add *tess/testthat/Rplots.pdf* to *.gitignore*. 

# EWCE 1.7.1 

## New features

* Use `rworkflows` GHA.
    - Add `rworkflows::use_badges` to *README.Rmd*.
    - Remove *Dockerfile* (no longer necessary).
    - Make all 3 platforms (Linux, Mac, Windows) use Bioc dev, 
        as `ewceData (>=1.7.1)` is now required, due to a fix made only 
        in the development version of `rtracklayer`.
* Remove `cowplot` dependency.
* Replace all `%>%` with `|>`

## Bug fixes

* `calc_quantiles`:
    - This function was only used in `filter_variance_quantiles` 
    - Compare `stats::ecdf` vs. `dplyr::ntile` methods.
    - Remove from `EWCE` as it's not longer used anywhere.
* `bin_columns_into_quantiles`:
    - Rename arg `matrixIn` --> `vec` to reflect what the function actually does.
* `filter_variance_quantiles`:
    - Change to use `bin_columns_into_quantiles` instead of `calc_quantiles` 
        to be consistent with how quantiles are handled in the rest of `EWCE`.
    - Updated tests in *test-get_celltype_table.r* to reflect that the 
        number of genes filtered is unaffected by the normalization procedure 
        (when quantiles are computed with `stats::quantile`).
* `ewce_plot`:
    - Celltypes were producing NAs because the names in the results
        were not always standardized in the same way as the CTD. 
        Now this is done internally.
    - Celltypes were not ordered factors, 
        meaning the dendrogram didn't line up correctly.
    - Switched from `gridArrange`/`cowplot` to `patchwork`.
    - Added dedicated unit tests file: *test-ewce_plot.r*

# EWCE 1.7.1

## New features

* Offline runs enabled with functions using reference datasets 
(from `ewceData`). These functions have the parameter `localhub` added to 
control this.

# EWCE 1.5.8

## Bug fixes

* GHA fix.

# EWCE 1.5.7

## Bug fixes

* orthogene dependency has been replacing user entered background gene list with
one generated from all known genes when species across gene lists and reference 
dataset are the same. This has now been fixed.

# EWCE 1.5.5

## Bug fixes

* `standardise_ctd`:
    - Always force "specificity_quantiles" to be one of the matrices 
    in each level. 

# EWCE 1.5.4

## New features

* `filter_ctd_genes`
    - Now exported.
    - Can handle standardized CTD format.
* `get_ctd_matrix_names`: New function to get a list of all data matrices in CTD.
    
## Bug fixes

* `check_ewce_genelist_inputs`: 
    - User reported potential bug in code:https://github.com/NathanSkene/EWCE/issues/71 
    - Fixed by removing conditional and instead always filtering out genes not present in CTD/SCT. 
* `standardise_ctd`:
    - Add `check_species()`
    - Ensure all matrices become sparse when `as_sparse=TRUE`.
    - Generalize to matrices of any name. 
* `fix_celltype_names`:
    - Ensure all celltype names are unique after standardization. 

# EWCE 1.5.3

## New features

* `genelistSpecies` now passed to `prepare_genesize_control_network` in 
`bootstrap_enrichment_test` meaning gene list species will be inferred from user
input.

# EWCE 1.5.2

## New features

* `drop_uninformative_genes`:
    - Expose new args: `dge_method`, `dge_test`, `min_variance_decile`
* `merged_ctd`: Actually merge the CTDs into one when `as_SCE=FALSE`.

## Bug fixes

* Remove hard-coded file path separators
(e.g. `sprintf("%s/MRK_List2.rpt", tempdir())`) to be more compatible with Windows.

# EWCE 1.5.1

## New features

* Made substantial updates to `orthogene`, so going through and making sure everything still works / is able to take advantage of new features (e.g. separation of `non121_strategy` and `agg_func` args, many:many mapping):  
    - `filter_nonorthologs`: Pass up args from `orthogene::convert_orthologs`.
    - `generate_celltype_data`: @inheritDotParams
* Update GHA. 
* Bump to R (>= 4.2) now that we're developing on Bioc 3.16. 

## Bug fixes

* Avoid downloading large "MRK_List2.rpt" file any more than 
is necessary for testing. 

# EWCE 1.3.3

## New features

* `method` argument from `orthogene::create_background` and `orthogene::convert_orthologs` is now passed up as an argument to `EWCE` functions to give users more control. "homologene" chosen as default for all functions. "homologene" has fewer species than "orthogene" but doesnt need to import data from the web. It also has more 1:1 mouse:human orthologs. 
* Include notes on mismatches between GitHub documentation and current Bioc release version. 
* Allow `bin_specificity_into_quantiles` to set specificity matrix
name produced. 
* Merge GHA workflow yamls into one. 

## Bug fixes

* Add `try({})` and `error=TRUE` to avoid *"polygon edge not found"* error in vignettes.  


# EWCE 1.3.1

## New features  

* Major changes: Pull Request from *bschilder_dev* branch. 
* All functions can now use lists and CellTypeDatasets (CTD) from any species 
and convert them to a common species (human by default) via `orthogene`.  
* Automated CTD standardisation via `standardise_ctd`.  
* Can handle (sparse) matrices.  
* Can create CTD from very large datasets using `DelayedArray` object class.  
* All functions automatically create appropriate gene backgrounds given species.  
* More modular, simplified vignettes.  
* Additional gene pre-filtering options (DESeq2, MAST, variance quantiles).  
* New/improved plotting functions (e.g. `plot_ctd`).  
* Added example bootstrapping enrichment results as *extdata* to 
speed up examples (documented in *data.R*). 
Accessed via `EWCE::example_bootstrap_results()`.
* Replaced GHA workflow with check-bioc to automatically: run R-CMD checks, 
run `BiocCheck`, and rebuild/deploy pkgdown site. 
* Parallelised functions: 
    + `drop_uninformative_genes`  
    + `generate_celltype-data`  
    + `bootstrap_enrichment_test`   
* Added tests (multiple functions tests per file to reduce number
of times `ewceData` files have to be downloaded):
    + `test-DelayedArray` 
    + `test-merge_sce`
    + `test-get_celltype_table`
    + `test-list_species` 
    + `test-run_DGE` 
    + `test-check_percent_hits` 
* Added function `is_32bit()` to all tests to ensure they don't 
get run twice on Windows OS.  
* Added GitHub Actions workflows:
    + `check-bioc-docker.yml`: Runs CRAN/Bioc checks, rebuilds and pushes `pkgdown` website, runs and uploads test coverage report, 
    + `dockerhub.yml`: Builds Bioconductor Docker container with `EWCE` installed, runs CRAN checks and (if checks are successful) pushes container to [*neurogenomicslab* DockerHub](https://hub.docker.com/repository/docker/neurogenomicslab/ewce).  
* Removed `docs` folder, as the documentation website comes from the
*gh-pages* branch now, and is automatically built by GHA workflow 
after each push to *main* branch. 
* Added new exported function `fix_celltype_names` to help with standardising 
celltype names in alignment with `standardise_ctd`. 
* `generate_bootstrap_plots_for_transcriptome`: Now supports any species 
(not just mouse or human). 
    + Converts CTD and DGE table (`tt`) into `output_species` gene symbols. 
    + Automatically generates appropriate gene background.  
    + Faster due to now having the option to only generate certain plot types.  
* Provide precomputed results from `ewce_expression_data` via new `example_transcriptome_results` function.  
* Reduced build runtime and oversized vignettes by not evaluating 
certain code chunks.  
    - Prevent *extended* vignette from running entirely. 
* `@return` documentation for internal functions.  
* Added more installation checks to GHA. 
* Fixed inconsistent naming of unit test files: `test_` ==> `test-` 
* Removed DGE args in `drop_uninformative_genes` for now until we run
benchmarking to see how each affects the `EWCE` results.  
* Make `bootstrap_plots` function internal. 
* Add report on how `orthogene` improve within- and across-species
gene mappings in *extended* vignette. 
* Record extra info in `standardise_ctd` output:
    - "species": both `input_species` and `output_species` 
    - "versions": of `EWCE`, `orthogene`, and `homologene`  


# EWCE 1.0.0

## New features  

* EWCE v1.0 on Bioconductor replaces the defunct [EWCE v1.3.0](https://bioconductor.riken.jp/packages/3.5/bioc/html/EWCE.html)
available on Bioconductor v3.5.
* EWCE has been rendered scalable to the analysis of large datasets
* `drop_uninformative_genes()` has been expanded to allow the utilisation of differential expression approaches 
* EWCE can now handle SingleCellExperiment (SCE) objects or other Ranged SummarizedExperiment (SE) data types and as input as well as the original format, described as a single cell transcriptome (SCT) object.

Deprecated & Defunct  

* The following functions have been renamed to use underscore in compliance with Bioconductor nomenclature:
    + `check.ewce.genelist.inputs`  
    + `cell.list.dist`  
    + `bootstrap.enrichment.test`  
    + `bin.specificity.into.quantiles`  
    + `bin.columns.into.quantiles`  
    + `add.res.to.merging.list`  
    + `prepare.genesize.control.network`  
    + `prep.dendro`  
    + `get.celltype.table`  
    + `calculate.specificity.for.level`  
    + `calculate.meanexp.for.level`  
    + `generate.celltype.data`  
    + `generate.bootstrap.plots`  
    + `generate.bootstrap.plots.for.transcriptome`  
    + `fix.bad.mgi.symbols`  
    + `fix.bad.hgnc.symbols`  
    + `filter.genes.without.1to1.homolog`  
    + `ewce.plot`  
    + `cells.in.ctd`  
    + `drop.uninformative.genes`  
