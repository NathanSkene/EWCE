`EWCE`: `E`xpression `W`eighted `C`elltype `E`nrichment
================
<img height='200' src='https://github.com/NathanSkene/EWCE/blob/bschilder_dev/inst/hex/hex.png?raw=true'><br><br>
[![](https://img.shields.io/badge/devel%20version-1.3.3-black.svg)](https://github.com/NathanSkene/EWCE)
[![](https://img.shields.io/badge/release%20version-1.2.0-green.svg)](https://www.bioconductor.org/packages/EWCE)
[![R build
status](https://github.com/NathanSkene/EWCE/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/NathanSkene/EWCE/actions)
[![](https://img.shields.io/github/last-commit/NathanSkene/EWCE.svg)](https://github.com/NathanSkene/EWCE/commits/master)
[![](https://codecov.io/gh/NathanSkene/EWCE/branch/master/graph/badge.svg)](https://codecov.io/gh/NathanSkene/EWCE)
[![](https://img.shields.io/badge/download-2422/total-green.svg)](https://bioconductor.org/packages/stats/bioc/EWCE)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
<h4>
Authors: <i>Alan Murphy, Brian Schilder, Nathan Skene</i>
</h4>
<h4>
README updated: <i>Feb-09-2022</i>
</h4>

<!-- To modify Package/Title/Description/Authors fields, edit the DESCRIPTION file -->

## Introduction

The *EWCE* R package is designed to facilitate expression weighted cell
type enrichment analysis as described in our *Frontiers in Neuroscience*
paper.<sup>1</sup> *EWCE* can be applied to any gene list.

Using *EWCE* essentially involves two steps:

1.  Prepare a single-cell reference; i.e. CellTypeDataset (CTD).
    Alternatively, you can use one of the pre-generated CTDs we provide
    via the package `ewceData` (which comes with *EWCE*).  
2.  Run cell type enrichment on a user-provided gene list.

## Installation

*EWCE* requires [`R>=4.1`](https://www.r-project.org/) and
`Bioconductor>=3.14`. To install *EWCE* on Bioconductor run:

``` r
if (!require("BiocManager")){install.packages("BiocManager")}

BiocManager::install("EWCE") 
```

## Documentation

### [Website](https://nathanskene.github.io/EWCE/)

**NOTE**: This documentation is for the development version of `EWCE`.
See
[Bioconductor](https://bioconductor.org/packages/release/bioc/html/EWCE.html)
for documentation on the current release version.

### [Getting started](https://nathanskene.github.io/EWCE/articles/EWCE)

Includes:

-   A minimal example to get started with running *EWCE*.
-   How to install and use the dedicated *EWCE* Docker container usage.
    [Docker](https://www.docker.com/) containers with the latest version
    of `EWCE` are regularly pushed to
    [Dockerhub](https://hub.docker.com/repository/docker/neurogenomicslab/ewce).

### [Extended examples](https://nathanskene.github.io/EWCE/articles/extended.html)

Additional tutorials of various *EWCE* features, including how to:

-   Run cell-type enrichment tests
-   Create a CellTypeDataset
-   Merge two single-cell datasets
-   Run conditional cell-type enrichment tests
-   Apply to transcriptomic data

## Updates

Major upgrades to *EWCE* were made in version 1.3.1. Please see the
[NEWS page](https://nathanskene.github.io/EWCE/news/index.html) for more
details.

## Troubleshooting

If you have any problems, please do submit an [Issue here on
GitHub](https://github.com/nathanskene/EWCE/issues) with a reproducible
example.

## Citation

If you use `EWCE`, please cite:

<!-- Modify this my editing the file: inst/CITATION  -->

> Nathan G. Skene, Seth G. N. Grant (2016) Identification of Vulnerable
> Cell Types in Major Brain Disorders Using Single Cell Transcriptomes
> and Expression Weighted Cell Type Enrichment, *Frontiers in
> Neuroscience*; 10, <https://doi.org/10.3389/fnins.2016.00016>

If you use the cortex/hippocampus single-cell data associated
*EWCE*/*ewceData* this package then please cite the following:

> [Zeisel, et al. Cell types in the mouse cortex and hippocampus
> revealed by single-cell RNA-seq. Science,
> 2015.](https://doi.org/10.1126/science.aaa1934)

<hr>

## Contact

### [Neurogenomics Lab](https://www.neurogenomics.co.uk/)

UK Dementia Research Institute  
Department of Brain Sciences  
Faculty of Medicine  
Imperial College London  
[GitHub](https://github.com/neurogenomics)  
[DockerHub](https://hub.docker.com/orgs/neurogenomicslab)

## References

<div id="refs" class="references csl-bib-body" line-spacing="2">

<div id="ref-skene_2016" class="csl-entry">

<span class="csl-left-margin">1. </span><span
class="csl-right-inline">Skene, N. & Grant, S. Identification of
vulnerable cell types in major brain disorders using single cell
transcriptomes and expression weighted cell type enrichment. *Frontiers
in Neuroscience* (2016).
doi:[10.3389/fnins.2016.00016](https://doi.org/10.3389/fnins.2016.00016)</span>

</div>

</div>
