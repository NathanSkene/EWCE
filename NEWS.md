# EWCE 2.0.0

New features  

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
run BiocCheck, and rebuild/deploy pkgdown site.


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
