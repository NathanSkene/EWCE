# EWCE 1.3.3

New features

* `method` argument from `orthogene::create_background` and `orthogene::convert_orthologs` is now passed up as an argument to `EWCE` functions to give users more control. "homologene" chosen as default for all functions. "homologene" has fewer species than "orthogene" but doesnt need to import data from the web. It also has more 1:1 mouse:human orthologs. 
* Include notes on mismatches between GitHub documentation and current Bioc release version. 

# EWCE 1.3.1

New features  

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

New Features  

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
