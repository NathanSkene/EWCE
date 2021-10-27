`EWCE`: `E`xpression `W`eighted `C`elltype `E`nrichment
================
<h4>
Alan Murphy, Brian Schilder, and Nathan Skene
</h4>
<h4>
Oct-27-2021
</h4>

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-2.0.0-black.svg)](https://github.com/NathanSkene/EWCE)
[![R build
status](https://github.com/NathanSkene/EWCE/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/NathanSkene/EWCE/actions)
[![R build
status](https://github.com/NathanSkene/EWCE/workflows/DockerHub/badge.svg)](https://github.com/NathanSkene/EWCE/actions)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/EWCE.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/EWCE)
[![platforms](http://www.bioconductor.org/images/shields/availability/all.svg)](https://bioconductor.org/packages/devel/bioc/html/EWCE.html#archives)
[![](https://img.shields.io/badge/doi-10.18129/B9.bioc.EWCE%20-green.svg)](https://doi.org/10.18129/B9.bioc.EWCE)
[![](https://img.shields.io/github/last-commit/NathanSkene/EWCE.svg)](https://github.com/NathanSkene/EWCE/commits/master)
[![](https://codecov.io/gh/NathanSkene/EWCE/branch/master/graph/badge.svg)](https://codecov.io/gh/NathanSkene/EWCE)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![](https://img.shields.io/badge/download-2155/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/EWCE)
<!-- badges: end -->

<img height="200" src="https://github.com/bschilder/EWCE/raw/DelayedArray/inst/hex/EWCE.png">

# Introduction

The *EWCE* R package is designed to facilitate expression weighted cell
type enrichment analysis as described in our *Frontiers in Neuroscience*
paper.<sup>1</sup> *EWCE* can be applied to any gene list.

Using *EWCE* essentially involves two steps:

1.  Prepare a single-cell reference; i.e. CellTypeDataset (CTD).
    Alternatively, you can use one of the pre-generated CTDs we provide
    via the package `ewceData` (which comes with *EWCE*).  
2.  Run cell type enrichment on a user-provided gene list.

# Installation

*EWCE* requires [`R>=4.1`](https://www.r-project.org/) and
`Bioconductor>=1.4`. To install *EWCE* on Bioconductor run:

``` r
if (!require("BiocManager")){install.packages("BiocManager")}

BiocManager::install(version = "devel")
BiocManager::install("EWCE") 
```

# [Documentation site](https://nathanskene.github.io/EWCE)

# Vignettes

## [EWCE: Getting started](https://nathanskene.github.io/EWCE/articles/EWCE.html)

A minimal example to get started with running *EWCE*.

## [EWCE: Extended examples](https://nathanskene.github.io/EWCE/articles/extended.html)

Example *EWCE* enrichment tests with in-depth explanations of each step.

## [EWCE: Creating CellTypeDatasets](https://nathanskene.github.io/EWCE/articles/create_CTD.html)

Instructions on how to create new make new CellTypeDataset references to
use with *EWCE*.

## [EWCE: Conditional examples](https://nathanskene.github.io/EWCE/articles/conditional.html)

Examples and explanations of conditional cell type enrichment tests
(e.g. controlling for a dominant cell type signal) .

## [EWCE: Transcriptome examples](https://nathanskene.github.io/EWCE/articles/transcriptomes.html)

Additional applications of *EWCE* to transcriptomic studies.

# Docker container

[Docker](https://www.docker.com/) containers with the latest version of
`EWCE` are regularly pushed to [Dockerhub](https://hub.docker.com/). If
you already have Docker installed, you can load up a working copy using
the following commands.

Note, that you will need to replace the initial directory path with a
location on your computer that you wish to be able to access from within
the docker image.

    docker pull nathanskene/ewce
    docker run --name=ewce -e PASSWORD=ewcedocker -p 8790:8790 -d -v /User/$USER:/var/ewce nathanskene/ewce:latest
    docker exec -ti ewce R

# Troubleshooting

If you have any problems, please do submit an [Issue here on
GitHub](https://github.com/nathanskene/EWCE/issues) with a reproducible
example.

# Citation

If you use EWCE, please cite:

> [Skene, et al. Identification of Vulnerable Cell Types in Major Brain
> Disorders Using Single Cell Transcriptomes and Expression Weighted
> Cell Type Enrichment. Front. Neurosci,
> 2016.](https://www.frontiersin.org/articles/10.3389/fnins.2016.00016/full)

If you use the cortex/hippocampus single-cell data associated
`EWCE`/`ewceData` this package then please cite the following:

> [Zeisel, et al. Cell types in the mouse cortex and hippocampus
> revealed by single-cell RNA-seq. Science,
> 2015.](http://www.sciencemag.org/content/early/2015/02/18/science.aaa1934.abstract)

# Session Info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.2 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] BiocManager_1.30.16 compiler_4.1.1      pillar_1.6.4       
    ##  [4] RColorBrewer_1.1-2  yulab.utils_0.0.4   tools_4.1.1        
    ##  [7] digest_0.6.28       jsonlite_1.7.2      evaluate_0.14      
    ## [10] lifecycle_1.0.1     tibble_3.1.5        gtable_0.3.0       
    ## [13] pkgconfig_2.0.3     rlang_0.4.12        DBI_1.1.1          
    ## [16] rvcheck_0.2.1       yaml_2.2.1          xfun_0.27          
    ## [19] fastmap_1.1.0       stringr_1.4.0       dplyr_1.0.7        
    ## [22] knitr_1.36          desc_1.4.0          generics_0.1.1     
    ## [25] vctrs_0.3.8         dlstats_0.1.4       rprojroot_2.0.2    
    ## [28] grid_4.1.1          tidyselect_1.1.1    glue_1.4.2         
    ## [31] R6_2.5.1            fansi_0.5.0         rmarkdown_2.11     
    ## [34] ggplot2_3.3.5       purrr_0.3.4         badger_0.1.0       
    ## [37] magrittr_2.0.1      scales_1.1.1        ellipsis_0.3.2     
    ## [40] htmltools_0.5.2     assertthat_0.2.1    colorspace_2.0-2   
    ## [43] utf8_1.2.2          stringi_1.7.5       munsell_0.5.0      
    ## [46] crayon_1.4.1

</details>

# References

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
