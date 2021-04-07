## CHANGES IN VERSION 1.0.0

### New Features

*   EWCE v1.0 on Bioconductor replaces the defunct [EWCE v1.3.0](https://bioconductor.riken.jp/packages/3.5/bioc/html/EWCE.html) available on Bioconductor v3.5.
*   EWCE has been rendered scalable to the analysis of large datasets
*   `drop_uninformative_genes()` has been expanded to allow the utilisation of differential expression approaches 
*   EWCE can now handle SingleCellExperiment (SCE) objects or other Ranged SummarizedExperiment (SE) data types and as input as well as the original format, described as a single cell transcriptome (SCT) object.

### Deprecated & Defunct

*   The following functions have been renamed to use underscore in compliance with Bioconductor nomenclature:
    `check.ewce.genelist.inputs`  
    `cell.list.dist`  
    `bootstrap.enrichment.test`  
    `bin.specificity.into.quantiles`  
    `bin.columns.into.quantiles`  
    `add.res.to.merging.list`  
    `prepare.genesize.control.network`  
    `prep.dendro`  
    `get.celltype.table`  
    `calculate.specificity.for.level`  
    `calculate.meanexp.for.level`  
    `generate.celltype.data`  
    `generate.bootstrap.plots`  
    `generate.bootstrap.plots.for.transcriptome`  
    `fix.bad.mgi.symbols`  
    `fix.bad.hgnc.symbols`  
    `filter.genes.without.1to1.homolog`  
    `ewce.plot`  
    `cells.in.ctd`  
    `drop.uninformative.genes`  
