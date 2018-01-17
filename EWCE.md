Expression Weighted Celltype Enrichment with *EWCE*
================
Nathan Skene
2018-01-17

<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>
Citation
========

If you use *[EWCE](http://bioconductor.org/packages/EWCE)*<sup>1</sup> in published research, please cite N. Skene et al (2016).

    N.G. SKENE, S.G.N. GRANT
    "Identification of Vulnerable Cell Types in Major Brain Disorders Using Single Cell Transcriptomes and Expression Weighted Cell Type Enrichment."
    Frontiers in Neuroscience, 27th Jan 2016

URL: <http://journal.frontiersin.org/article/10.3389/fnins.2016.00016/abstract>

Introduction
============

The *EWCE* package is designed to facilitate expression weighted celltype enrichment analysis as described in our Frontiers in Neuroscience paper<sup>1</sup>.

The package was originally designed to work with the single cell cortical transcriptome data from the Linnarsson lab<sup>2</sup> which is available at <http://linnarssonlab.org/cortex/>. Using this package it is possible to read in any single cell transcriptome data, provided that you have a cell by gene expression matrix (with each cell as a seperate column) and a seperate annotation dataframe, with a row for each cell.

The *EWCE* process involves testing for whether the genes in a target list have higher levels of expression in a given cell type than can reasonably be expected by chance. The probability distribution for this is estimated by randomly generating gene lists of equal length from a set of background genes.

The *EWCE* method can be applied to any gene list. In the paper we reported it's application to genetic and transcriptomic datasets, and in this vignette we detail how this can be done.

Note that throughout this vignette we use the terms 'cell type' and 'sub-cell type' to refer to two levels of annotation of what a cell type is. This is described in further detail in our paper<sup>1</sup>, but relates to the two levels of annotation provided in the Linnarsson dataset<sup>2</sup>. In this dataset a cell is described as having a cell type (i.e. 'Interneuron') and subcell type (i.e. 'Int11' a.k.a Interneuron type 11).

Overview
========

The process for using *EWCE* essentially involves three steps.

First, one needs to load the relevant single cell transcriptome dataset. Single cell transcriptome data is read in from a text file using the `read_celltype_data`.

The user then obtains a gene set and a suitable background gene set. As the choice of gene sets is up to the user we do not provide functions for doing this. Appropriate choice of background set is discussed in the associated publication.

Bootstrapping is then performed using the `bootstrap.enrichment.test` function.

Installing EWCE
===============

The *EWCE* package is available from the Bioconductor repository at <http://www.bioconductor.org> To be able to install the package one needs first to install R and the core Bioconductor packages. If you have already installed Bioconductor packages on your system then you can skip the two lines below.

    source("http://bioconductor.org/biocLite.R")
    biocLite()

Once the core Bioconductor packages are installed, we can install the *EWCE* package by

    source("http://bioconductor.org/biocLite.R")
    biocLite("EWCE")

You can then load the package:

``` r
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)
```

Loading single cell transcriptome data
======================================

Loading datasets
----------------

The first step for all analyses is to load the single cell transcriptome (SCT) data. For the purposes of this example we will use the dataset described in *"Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq", Science, 2015*. The following code downloads this file, passes it to the load.linnarsson.sct.data() function which extracts the expression and annotation data and returns these as a list.

Important note: you do NOT have to format your input single cell data like the Linnarsson data. Just read it into R such that you have an expression matrix and an annotation data frame. The three columns that you must have in the annotation data frame are "cell\_id", "level1class" and "level2class".

    download.file("goo.gl/r5Y24y",
        destfile="expression_mRNA_17-Aug-2014.txt") 

    path = "expression_mRNA_17-Aug-2014.txt"

    cortex_mrna  = load.linnarsson.sct.data(path)

To check the data we can quickly plot the distribution of expression of a given gene across all the cell types.

``` r
data("cortex_mrna")
gene="Necab1"
cellExpDist = data.frame(e=cortex_mrna$exp[gene,],l1=cortex_mrna$annot[colnames(cortex_mrna$exp),]$level1class)
ggplot(cellExpDist) + geom_boxplot(aes(x=l1,y=e)) + xlab("Cell type") + ylab("Unique Molecule Count") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-2-1.png)

It is very common for publically available transcriptome datasets to use incorrect gene symbols (often some gene names will have been mangled by opening in Excel). A function is provided to correct these where out of date aliases were used and give a warning if appears possible that excel has mangled some of the gene names. The first time you run it you will need to download the MRK\_List2 file from MGI which lists known synonyms for MGI gene symbols. We recommend running this on all input datasets.

    if(!file.exists("MRK_List2.rpt")){
        download.file("http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt", destfile="MRK_List2.rpt")
    }
    cortex_mrna$exp = fix.bad.mgi.symbols(cortex_mrna$exp,mrk_file_path="MRK_List2.rpt")

Calculate specificity matrices
------------------------------

The next steps are as follows: 1) Drop genes which do not show significant evidence of varying between level 2 celltypes (based on ANOVA) 2) Calculate cell type averages and specificity for each gene 3) Drop all genes which do not have 1:1 mouse:human orthologs

The last step is only neccesary if you plan to compare to human data, i.e. genesets resulting from human genetics.

Rather than returning the data directly, the functions save the calculated data to a file and the filename is returned.

``` r
# Generate celltype data for just the cortex/hippocampus data
exp_CortexOnly_DROPPED = drop.uninformative.genes(exp=cortex_mrna$exp,level2annot = cortex_mrna$annot$level2class)
fNames_CortexOnly = generate.celltype.data(exp=exp_CortexOnly_DROPPED,level1class=cortex_mrna$annot$level1class,level2class=cortex_mrna$annot$level2class,groupName="kiCortexOnly",thresh=0,trim=0)
print(fNames_CortexOnly)
```

    ## [1] "CellTypeData_kiCortexOnly.rda"

``` r
fNames_CortexOnly = filter.genes.without.1to1.homolog(fNames_CortexOnly)
print(fNames_CortexOnly)
```

    ## [1] "CellTypeData_kiCortexOnly.rda"         
    ## [2] "CellTypeData_kiCortexOnly_1to1only.rda"

``` r
load(fNames_CortexOnly[1])
```

Merging two single cell datasets
--------------------------------

Often it is useful to merge two single cell datasets. For instance, there are seperate files for the cortex and hypothalamus datasets generated by the Karolinska. This dataset is first downloaded from GEO, unzipped and read into R. Because the file is in xlsx format you may need to install the readxl package. So first we just read in and prepare the expression matrix and annotation data frame:

    # Download the hypothalamus data and unzip
    if(!file.exists("GSE74672_expressed_mols_with_classes.xlsx")){
        download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74672/suppl/GSE74672_expressed_mols_with_classes.xlsx.gz", destfile="GSE74672_expressed_mols_with_classes.xlsx.gz")
        system("gunzip GSE74672_expressed_mols_with_classes.xlsx.gz")
    }

    # Read in the hypothalamus data
    hypo_dat = read_excel("GSE74672_expressed_mols_with_classes.xlsx")

    # Extract the expression data, gene symbols and annotation data
    exp = data.matrix(hypo_dat[12:dim(hypo_dat)[1],2:dim(hypo_dat)[2]])
    rownames(exp) = data.frame(hypo_dat[12:dim(hypo_dat)[1],1])[,1]
    level1class = data.frame(level1class=t(hypo_dat[1,2:dim(hypo_dat)[2]]),stringsAsFactors = FALSE)[,1]
    level2class = data.frame(leve2class=t(hypo_dat[2,2:dim(hypo_dat)[2]]),stringsAsFactors = FALSE)[,1]
    cell_id     = colnames(hypo_dat)[2:dim(hypo_dat)[2]]
    hypo_annot  = data.frame(cell_id=cell_id,level1class=level1class,level2class=level2class,stringsAsFactors = FALSE)

    # Drop the glia and unclassified cells (which don't have level 2  annotations)
    hypo_annot  = hypo_annot[!is.na(hypo_annot$level2class) & !hypo_annot$level2class=="uc",]
    hypo_exp    = exp[,hypo_annot$cell_id]

    # Make the celltype names more aesthetically pleasing
    hypo_annot$level2class=gsub(",",";",hypo_annot$level2class)
    hypo_annot$level1class[grep("Oxt;|^Avp",hypo_annot$level2class)] = "Oxytocin / Vasopressin Expressing Neurons"
    hypo_annot$level1class[grep("^Th;|^Dopamine",hypo_annot$level2class)] = "Hypothalamic Dopaminergic Neurons"
    hypo_annot$level1class[grepl("^Vglut2|^Trh|^Qrfp|^Hcrt|^Pmch|^Adcyap1|^Npvf|^Ghrh|^Hmit|^Nms|^Vip;|^Per2|Tnr$|^Gad-low;Gnrh",hypo_annot$level2class) & grepl("neurons",hypo_annot$level1class)] = "Hypothalamic Glutamatergic Neurons"
    hypo_annot$level1class[grepl("GABA|^Sst|^Crh|^Npy|^Pomc|^Galanin|^Otof|Pnoc$|^Calcr-high",hypo_annot$level2class) & grepl("^neurons$",hypo_annot$level1class)] = "Hypothalamic GABAergic Neurons"
    hypo_annot$level2class[hypo_annot$level2class!=""] = sprintf("Hypothalamic %s Neuron",hypo_annot$level2class[hypo_annot$level2class!=""])

    # Fix bad MGI symbols
    hypo_exp_CORRECTED = fix.bad.mgi.symbols(hypo_exp)

Now that the hypothalamus data is prepared, we merge it with the cortex dataset then calculate specificity:

    # Merge the datasets
    merged_KI = merge_two_expfiles(exp1=hypo_exp_CORRECTED,  exp2=cortex_mrna$exp,
                                         annot1=hypo_annot,        annot2=cortex_mrna$annot,
                                         name1="Hypothalamus (KI)", name2="Cortex/Hippo (KI)")

    # Drop genes which don't vary significantly between cell types
    exp_merged_DROPPED = drop.uninformative.genes(exp=merged_KI$exp, level2annot = merged_KI$annot$level2class)

    # Calculate specificity data
    fNames_MergedKI = generate.celltype.data(exp=exp_merged_DROPPED,merged_KI$annot$level1class,merged_KI$annot$level2class,"MergedKI",thresh=0,trim=0)
    fNames_MergedKI = filter.genes.without.1to1.homolog(fNames_MergedKI)
    load(fNames_MergedKI[2])

Understanding specificity matrices
----------------------------------

While not required for further analyses it helps to understand what the outputs of this function are.

Note firstly, that it is a list such that ctd\[\[1\]\] contains data relating to level 1 annotations and ctd\[\[2\]\] relates to level 2 annotations.

Using the ggplot2 package to visualise the data, let us examine the expression of a few genes. If you have not already done so you will need to first install the ggplot2 package with `install.packages("ggplot2")`.

For this example we use a subset of the genes from the merged dataset generated above, which is accessed using `data(ctd)`. We recommend that you use the code above to regenerate this though and drop the `data` command from the below section.

``` r
data("ctd")
set.seed(1234)
library(reshape2)
genes = c("Apoe","Gfap","Gapdh")
exp = melt(cbind(ctd[[1]]$mean_exp[genes,],genes),id.vars="genes")
colnames(exp) = c("Gene","Cell","AvgExp")
ggplot(exp)+geom_bar(aes(x=Cell,y=AvgExp),stat="identity")+facet_grid(Gene~.)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-4-1.png)

This graph shows the average expression of three genes: *Apoe, Gfap* and *Gapdh*. While there are substantial differences in which cell types express these genes, the dominant effect seen here is the overall expression level of the data. For the purposes of this analysis though, we are not interested in overall expression level and only wish to know about the proportion of a genes expression which is found in a particular celltype. We can study this instead using the following code which examines the data frame *ctd\[\[1\]\]$specificity*:

``` r
exp = melt(cbind(data.frame(ctd[[1]]$specificity[genes,]),genes),id.vars="genes")
colnames(exp) = c("Gene","Cell","Expression")
ggplot(exp)+geom_bar(aes(x=Cell,y=Expression),stat="identity")+facet_grid(Gene~.)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png)

We can now see in this graph that *Gfap* is the most specific to a cell type (Type 1 Astrocytes) of either of those three genes, with over 60% of it's expression found in that cell type.

It can also be seen that the majority of expression of Gapdh is in neurons but because their are a greater number of neuronal subtypes, the total expression proportion appears lower. We can examine expression across level 2 celltype level annotations by looking at *ctd\[\[2\]\]$specificity*:

``` r
exp = melt(cbind(data.frame(ctd[[2]]$specificity[genes,]),genes),id.vars="genes")
colnames(exp) = c("Gene","Cell","Specificity")
ggplot(exp)+geom_bar(aes(x=Cell,y=Specificity),stat="identity")+facet_grid(Gene~.)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-1.png)

Application to genetic data
===========================

Preparing gene lists
--------------------

For the first demonstration of EWCE we will test for whether genes that are genetically associated with Alzheimer's disease are enriched in any particular celltype. This gene list is stored within the package are we access it by first loading the package and then the dataset:

``` r
data("example_genelist")
print(example_genelist)
```

    ##  [1] "APOE"     "BIN1"     "CLU"      "ABCA7"    "CR1"      "PICALM"  
    ##  [7] "MS4A6A"   "CD33"     "MS4A4E"   "CD2AP"    "EOGA1"    "INPP5D"  
    ## [13] "MEF2C"    "HLA-DRB5" "ZCWPW1"   "NME8"     "PTK2B"    "CELF1"   
    ## [19] "SORL1"    "FERMT2"   "SLC24A4"  "CASS4"

All gene IDs are assumed by the package to be provided in gene symbol format (rather than Ensembl/Entrez). Symbols can be provided as either HGNC or MGI symbols, though the genelistSpecies argument will need to be set appropriately. Likewise, the single cell dataset can use either human or mouse gene symbols, but the sctSpecies argument must be set to either "human" or "mouse". The default species for both is mouse.

The example gene list here stores the human genes associated with disease, and hence are HGNC symbols.

The next step is to determine the most suitable background set. The experimental methods used to find these gene are all genome wide, so there is no restriction imposed as a result of that. Thus our initial background set is the set of all human genes. Not all human genes have mouse orthologs however, so we need to drop all genes from the target and background set which do not have mouse orthologs. To save repeatedly querying biomaRt we have a stored dataset containing all the human orthologs of MGI genes, `mouse_to_human_homologs`. We can use this to obtain the mouse orthologs of the target and background genes at the same time as we drop genes without orthologs:

``` r
data("mouse_to_human_homologs")
m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
mouse.hits = unique(m2h[m2h$HGNC.symbol %in% example_genelist,"MGI.symbol"])
#mouse.bg  = unique(setdiff(m2h$MGI.symbol,mouse.hits))
mouse.bg  = unique(m2h$MGI.symbol)
```

The target list is now converted to MGI symbols:

``` r
print(mouse.hits)
```

    ##  [1] "Apoe"    "Inpp5d"  "Cd2ap"   "Nme8"    "Cass4"   "Mef2c"   "Zcwpw1" 
    ##  [8] "Bin1"    "Clu"     "Celf1"   "Abca7"   "Slc24a4" "Ptk2b"   "Picalm" 
    ## [15] "Fermt2"  "Sorl1"

And we have 15604 genes in background set.

Setting analysis parameters
---------------------------

We now need to set the parameters for the analysis. For a publishable analysis we would want to generate over 10000 random lists and determine their expression levels, but for computational speed let us only use `reps=1000`. We want to analyse level 1 annotations so set level to 1.

``` r
reps=1000 # <- Use 1000 bootstrap lists so it runs quickly, for publishable analysis use >10000
level=1 # <- Use level 1 annotations (i.e. Interneurons)
```

Running EWCE analysis on genetic data
-------------------------------------

We have now loaded the SCT data, prepared the gene lists and set the parameters. We run the model as follows:

``` r
# Bootstrap significance testing, without controlling for transcript length and GC content
full_results = bootstrap.enrichment.test(sct_data=ctd,hits=mouse.hits,bg=mouse.bg,
                                reps=reps,annotLevel=level)
```

    ##  [1] "Apoe"    "Inpp5d"  "Cd2ap"   "Nme8"    "Cass4"   "Mef2c"   "Zcwpw1" 
    ##  [8] "Bin1"    "Clu"     "Celf1"   "Abca7"   "Slc24a4" "Ptk2b"   "Picalm" 
    ## [15] "Fermt2"  "Sorl1"  
    ## [1] "astrocytes_ependymal"
    ## [1] 0.054
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0.601
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 1
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0.014
    ## [1] "Fold enrichment: 2.02369844284128"
    ## [1] "Standard deviations from mean: 2.69149534380412"
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.524
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 0.647
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0.663
    ## [1] ""

The main table of results is stored in `full_results$results`. We can see the most significant results using:

``` r
print(full_results$results[order(full_results$results$p),3:5][1:6,])
```

    ##                          p fold_change sd_from_mean
    ## microglia            0.014   2.0236984    2.6914953
    ## astrocytes_ependymal 0.054   1.5960026    1.8111846
    ## oligodendrocytes     0.524   0.9465844   -0.2156873
    ## endothelial-mural    0.601   0.8461038   -0.4393394
    ## pyramidal CA1        0.647   0.9388163   -0.4008988
    ## pyramidal SS         0.663   0.9280526   -0.4991304

The results can be visualised using another function, which shows for each cell type, the number of standard deviations from the mean the level of expression was found to be in the target gene list, relative to the bootstrapped mean:

``` r
print(ewce.plot(full_results$results,mtc_method="BH"))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-13-1.png)

If you want to view the characteristics of enrichment for each gene within the list then the `generate.bootstrap.plots` function should be used. This saves the plots into the BootstrapPlots folder. This takes the results of a bootstrapping analysis so as to only generate plots for significant enrichments. The `listFileName` argument is used to give the generated graphs a particular file name.

``` r
generate.bootstrap.plots(sct_data=ctd,hits=mouse.hits,bg=mouse.bg,reps=100,annotLevel=1,full_results=full_results,listFileName="VignetteGraphs")
```

Running EWCE analysis on genetic data with controls for transcript length and GC-content
----------------------------------------------------------------------------------------

When analysing genes found through genetic association studies it is important to consider biases which might be introduced as a result of transcript length and GC-content. The package can control for these by selecting the bootstrap lists such that the *i<sup>th</sup>* gene in the random list has properties similar to the*i<sup>th</sup>* gene in the target list. To enable the algorithm to do this it needs to be passed the gene lists as HGNC symbols rather than MGI.

``` r
#human.hits = unique(m2h[m2h$HGNC.symbol %in% example_genelist,"HGNC.symbol"])
#human.bg = unique(setdiff(m2h$HGNC.symbol,human.hits))
human.hits = example_genelist
human.bg = unique(c(human.hits,m2h$HGNC.symbol))
```

The bootstrapping function then takes different arguments:

``` r
# Bootstrap significance testing controlling for transcript length and GC content
cont_results = bootstrap.enrichment.test(sct_data=ctd,hits=human.hits,
                    bg=human.bg,reps=reps,annotLevel=1,geneSizeControl=TRUE,genelistSpecies="human",sctSpecies="mouse")
```

    ## [1] "CONTROLLED BOOTSTRAPPING NETWORK GENERATED"
    ##  [1] "Abca7"   "Apoe"    "Bin1"    "Cass4"   "Cd2ap"   "Celf1"   "Clu"    
    ##  [8] "Fermt2"  "Inpp5d"  "Mef2c"   "Nme8"    "Picalm"  "Ptk2b"   "Slc24a4"
    ## [15] "Sorl1"   "Zcwpw1" 
    ## [1] "astrocytes_ependymal"
    ## [1] 0.013
    ## [1] "Fold enrichment: 2.06962086219617"
    ## [1] "Standard deviations from mean: 2.88687446000956"
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0.469
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 0.953
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0.005
    ## [1] "Fold enrichment: 2.45457643769155"
    ## [1] "Standard deviations from mean: 3.2031226729245"
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.29
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 0.147
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0.213
    ## [1] ""

We plot these results using `ewce.plot`:

``` r
print(ewce.plot(cont_results$results,mtc_method="BH"))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-17-1.png)

This shows that the controlled method generates enrichments that are generally comparable to the standard method.

Running EWCE analysis on cell type level annotations
----------------------------------------------------

Both the analyses shown above were run on level 1 annotations. It is possible to test on the level 2 cell type level annotations by changing one of the arguments.

``` r
# Bootstrap significance testing controlling for transcript length and GC content
cont_results = bootstrap.enrichment.test(sct_data=ctd,hits=human.hits,
                    bg=human.bg,reps=reps,annotLevel=2,geneSizeControl=TRUE,genelistSpecies="human",sctSpecies="mouse")
```

    ## [1] "CONTROLLED BOOTSTRAPPING NETWORK GENERATED"
    ##  [1] "Abca7"   "Apoe"    "Bin1"    "Cass4"   "Cd2ap"   "Celf1"   "Clu"    
    ##  [8] "Fermt2"  "Inpp5d"  "Mef2c"   "Nme8"    "Picalm"  "Ptk2b"   "Slc24a4"
    ## [15] "Sorl1"   "Zcwpw1" 
    ## [1] "(none)"
    ## [1] 0.114
    ## [1] ""
    ## [1] "Astro1"
    ## [1] 0.018
    ## [1] "Fold enrichment: 2.61473007787905"
    ## [1] "Standard deviations from mean: 2.88283187415692"
    ## [1] ""
    ## [1] "Astro2"
    ## [1] 0.038
    ## [1] "Fold enrichment: 2.45232633391661"
    ## [1] "Standard deviations from mean: 2.35168645455456"
    ## [1] ""
    ## [1] "CA1Pyr1"
    ## [1] 0.025
    ## [1] "Fold enrichment: 1.56521188647139"
    ## [1] "Standard deviations from mean: 2.23301073446451"
    ## [1] ""
    ## [1] "CA1Pyr2"
    ## [1] 0.107
    ## [1] ""
    ## [1] "CA1PyrInt"
    ## [1] 0.011
    ## [1] "Fold enrichment: 1.7503220199075"
    ## [1] "Standard deviations from mean: 2.83928093073704"
    ## [1] ""
    ## [1] "CA2Pyr2"
    ## [1] 0.273
    ## [1] ""
    ## [1] "Choroid"
    ## [1] 0.239
    ## [1] ""
    ## [1] "ClauPyr"
    ## [1] 0.562
    ## [1] ""
    ## [1] "Epend"
    ## [1] 0.044
    ## [1] "Fold enrichment: 3.02021582245981"
    ## [1] "Standard deviations from mean: 2.37400845572311"
    ## [1] ""
    ## [1] "Int1"
    ## [1] 0.821
    ## [1] ""
    ## [1] "Int10"
    ## [1] 0.423
    ## [1] ""
    ## [1] "Int11"
    ## [1] 0.886
    ## [1] ""
    ## [1] "Int12"
    ## [1] 0.659
    ## [1] ""
    ## [1] "Int13"
    ## [1] 0.947
    ## [1] ""
    ## [1] "Int14"
    ## [1] 0.569
    ## [1] ""
    ## [1] "Int15"
    ## [1] 0.682
    ## [1] ""
    ## [1] "Int16"
    ## [1] 0.821
    ## [1] ""
    ## [1] "Int2"
    ## [1] 0.87
    ## [1] ""
    ## [1] "Int3"
    ## [1] 0.444
    ## [1] ""
    ## [1] "Int4"
    ## [1] 0.864
    ## [1] ""
    ## [1] "Int5"
    ## [1] 0.86
    ## [1] ""
    ## [1] "Int6"
    ## [1] 0.753
    ## [1] ""
    ## [1] "Int7"
    ## [1] 0.911
    ## [1] ""
    ## [1] "Int8"
    ## [1] 0.6
    ## [1] ""
    ## [1] "Int9"
    ## [1] 0.948
    ## [1] ""
    ## [1] "Mgl1"
    ## [1] 0.057
    ## [1] ""
    ## [1] "Mgl2"
    ## [1] 0.008
    ## [1] "Fold enrichment: 4.69793330506261"
    ## [1] "Standard deviations from mean: 4.52115860930145"
    ## [1] ""
    ## [1] "Oligo1"
    ## [1] 0.4
    ## [1] ""
    ## [1] "Oligo2"
    ## [1] 0.155
    ## [1] ""
    ## [1] "Oligo3"
    ## [1] 0.229
    ## [1] ""
    ## [1] "Oligo4"
    ## [1] 0.147
    ## [1] ""
    ## [1] "Oligo5"
    ## [1] 0.183
    ## [1] ""
    ## [1] "Oligo6"
    ## [1] 0.182
    ## [1] ""
    ## [1] "Peric"
    ## [1] 0.157
    ## [1] ""
    ## [1] "Pvm1"
    ## [1] 0.022
    ## [1] "Fold enrichment: 3.21043003561403"
    ## [1] "Standard deviations from mean: 2.80724449283827"
    ## [1] ""
    ## [1] "Pvm2"
    ## [1] 0.082
    ## [1] ""
    ## [1] "S1PyrDL"
    ## [1] 0.284
    ## [1] ""
    ## [1] "S1PyrL23"
    ## [1] 0.024
    ## [1] "Fold enrichment: 1.77521499106251"
    ## [1] "Standard deviations from mean: 2.7605385741186"
    ## [1] ""
    ## [1] "S1PyrL4"
    ## [1] 0.095
    ## [1] ""
    ## [1] "S1PyrL5"
    ## [1] 0.158
    ## [1] ""
    ## [1] "S1PyrL5a"
    ## [1] 0.052
    ## [1] ""
    ## [1] "S1PyrL6"
    ## [1] 0.106
    ## [1] ""
    ## [1] "S1PyrL6b"
    ## [1] 0.164
    ## [1] ""
    ## [1] "SubPyr"
    ## [1] 0.567
    ## [1] ""
    ## [1] "Vend1"
    ## [1] 0.5
    ## [1] ""
    ## [1] "Vend2"
    ## [1] 0.328
    ## [1] ""
    ## [1] "Vsmc"
    ## [1] 0.411
    ## [1] ""

``` r
print(ewce.plot(cont_results$results,mtc_method="BH"))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-18-1.png)

With the subcell analysis each microglial subtype was enriched and correspondingly we see here that the microglial celltype is enriched.

Plotting results from multiple gene lists
-----------------------------------------

It is often useful to plot results from multiple gene list analyses together. The `ewce.plot` function allows multiple enrichment analyses to be performed together. To achieve this the results data frames are just appended onto each other, with an additional `list` column added detailing which analysis they relate to.

To demonstrate this we need to first generate a second analysis so let us sample thirty random genes, and run the bootstrapping analysis on it.

``` r
gene.list.2 = mouse.bg[1:30]
second_results = bootstrap.enrichment.test(sct_data=ctd,hits=gene.list.2,
                    bg=mouse.bg,reps=reps,annotLevel=1)
```

    ##  [1] "Rcl1"   "Myom1"  "Gpr12"  "Micu2"  "Pds5b"  "Jagn1"  "Rbm17" 
    ##  [8] "Trib2"  "Theg"   "Itpr1"  "Prmt5"  "Usp6nl" "Fgf9"   "Map1a" 
    ## [15] "Samd11" "Noc2l" 
    ## [1] "astrocytes_ependymal"
    ## [1] 0.93
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0.329
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 0.178
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0.868
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.702
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 0.387
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0.073
    ## [1] ""

``` r
full_res2 = data.frame(full_results$results,list="Alzheimers")
scnd_res2 = data.frame(second_results$results,list="Second")
merged_results = rbind(full_res2,scnd_res2)
```

``` r
print(ewce.plot(total_res=merged_results,mtc_method="BH"))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-20-1.png)

As expected, the second randomly generated gene list shows no significant enrichments.

Application to transcriptomic data
==================================

Analysing single transcriptome study
------------------------------------

For the prior analyses the gene lists were not associated with any numeric values or directionality. The methodology for extending this form of analysis to transcriptomic studies simply involves thresholding the most upregulated and downregulated genes.

To demonstrate this we have an example dataset `tt_alzh`. This data frame was generated using limma from a set of post-mortem tissue samples from brodmann area 46 which were described in a paper by the Haroutunian lab<sup>3</sup>.

The first step is to load the data, obtain the MGI ids, sort the rows by t-statistic and then select the most up/down-regulated genes. The package then has a function `ewce_expression_data` which thresholds and selects the gene sets, and calls the EWCE function. Below we show the function call using the default settings, but if desired different threshold values can be used, or alternative columns used to sort the table.

``` r
data(tt_alzh)
tt_results = ewce_expression_data(sct_data=ctd,tt=tt_alzh,annotLevel=1,ttSpecies="human",sctSpecies="mouse")
```

    ##   [1] "Myom1"         "Rbm17"         "Max"           "Mbd3"         
    ##   [5] "Wdr33"         "Acss1"         "Pcyt1b"        "Bcam"         
    ##   [9] "Prpf38b"       "Ppp2r5c"       "Ica1"          "Cobl"         
    ##  [13] "Marveld1"      "Gemin4"        "Zfp668"        "Casc3"        
    ##  [17] "Aff3"          "Tcf3"          "Taf1c"         "Gna12"        
    ##  [21] "Lamp2"         "Inpp5d"        "Phf14"         "Stk11"        
    ##  [25] "Plxnd1"        "Aebp1"         "Pgf"           "5430427O19Rik"
    ##  [29] "Gad1"          "App"           "Casp9"         "Mid1ip1"      
    ##  [33] "Cep70"         "Dnajc1"        "Cbx5"          "Llph"         
    ##  [37] "Rbm6"          "Tbc1d2b"       "Csk"           "Col6a1"       
    ##  [41] "1110034G24Rik" "Vars2"         "Lama2"         "Aup1"         
    ##  [45] "Ctbp2"         "Atg3"          "Aifm3"         "Gga1"         
    ##  [49] "Zcchc24"       "Prtg"          "Ralgps1"       "Ncan"         
    ##  [53] "Otud4"         "Phf21a"        "Adra2c"        "Iqcc"         
    ##  [57] "Ptbp1"         "Col27a1"       "Zfp526"        "Nrp2"         
    ##  [61] "Flnb"          "Dcn"           "Zfp36l1"       "Antxr1"       
    ##  [65] "Plod3"         "Ddit4l"        "Hnrnph3"       "Apba3"        
    ##  [69] "Taf5"          "Mmd2"          "Kif1c"         "Zswim1"       
    ##  [73] "Myh14"         "Wwox"          "Tmem47"        "Zfp36l2"      
    ##  [77] "Palld"         "Mfge8"         "Nfic"          "Ilf3"         
    ##  [81] "Armc3"         "Myh11"         "Ppp4r1"        "Plekha8"      
    ##  [85] "Gas2l1"        "Tal1"          "Atp8b2"        "Slc40a1"      
    ##  [89] "Enah"          "1110059E24Rik" "Nsd1"          "Rprm"         
    ##  [93] "Golim4"        "Paox"          "Numa1"         "Rabep1"       
    ##  [97] "Tcf12"         "Fastk"         "Cnot4"         "Cxcr4"        
    ## [101] "Rrad"          "Dmc1"          "Rassf9"        "Abca9"        
    ## [105] "Hbp1"          "Klk7"          "Nedd4"         "Adck1"        
    ## [109] "Smarcc1"       "Ezr"           "Gpr161"        "Rapgef3"      
    ## [113] "Snrnp48"       "Abi3bp"        "Nav1"          "Ssr1"         
    ## [117] "Pknox1"        "Qk"            "Osbp2"         "Col5a3"       
    ## [121] "Fgfr1"         "Repin1"        "Dlgap4"        "Atp5sl"       
    ## [125] "Rimkla"        "Pdcd2"         "Alcam"         "P2rx7"        
    ## [129] "Smyd4"         "Mtss1l"        "Sfxn5"         "Unc5b"        
    ## [133] "Ddc"           "Pcdh1"         "Igf2"          "Sec11a"       
    ## [137] "Osbpl2"        "3110079O15Rik" "Sall3"         "Git1"         
    ## [141] "Llgl2"         "Zfp385a"       "Kank2"         "Rps21"        
    ## [145] "Smg5"          "Gnao1"         "Csnk1d"        "Ubxn1"        
    ## [149] "Sos1"          "Ascl1"         "Map2k3"        "Nf2"          
    ## [153] "Fbxw4"         "Ocln"          "Zc3h14"        "S100a4"       
    ## [157] "Epb41l2"       "Rftn2"         "Zfp444"        "Rab28"        
    ## [161] "Notch2"        "Efcab2"        "Maf"           "Hes5"         
    ## [165] "Chdh"          "Mdm2"          "Arhgap22"      "Rnf11"        
    ## [169] "Crtap"         "Frzb"          "Gja1"          "Lmna"         
    ## [173] "Rheb"          "Tmem184b"      "Shroom1"       "Ahnak"        
    ## [177] "Itsn1"         "Jmjd6"         "Kcnj10"        "Bgn"          
    ## [181] "Atg12"         "Ngfr"          "Slc25a29"      "Pcmtd1"       
    ## [185] "Adam33"        "Rnd2"          "5031439G07Rik" "Mif4gd"       
    ## [189] "Atp6v0e"       "Hist3h2a"      "Tmcc2"         "Mrgprf"       
    ## [193] "Nrip2"         "Lef1"          "Sox12"         "Med25"        
    ## [197] "Jund"          "Acp6"          "Ak1"           "Cd79b"        
    ## [201] "Wwc1"          "Itga7"         "Fam208a"       "Map3k11"      
    ## [205] "Cspg4"         "Vat1"          "Zfp692"        "Selplg"       
    ## [209] "Dbi"           "Pou3f2"        "Doc2b"         "Itgb5"        
    ## [213] "Polr3g"        "Pacs1"        
    ## [1] "astrocytes_ependymal"
    ## [1] 0
    ## [1] "Fold enrichment: 1.38454551977423"
    ## [1] "Standard deviations from mean: 4.71834514988499"
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0
    ## [1] "Fold enrichment: 1.35882778152608"
    ## [1] "Standard deviations from mean: 4.72133515494952"
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 1
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0
    ## [1] "Fold enrichment: 1.32685008432062"
    ## [1] "Standard deviations from mean: 3.47950829305971"
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.24
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 1
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 1
    ## [1] ""
    ##   [1] "Prmt5"         "Atg7"          "Ivd"           "Ica1"         
    ##   [5] "Farsb"         "Aifm1"         "Pole3"         "Psmd11"       
    ##   [9] "Dpm2"          "Vrk3"          "Dnajc5"        "Gphn"         
    ##  [13] "Tubb5"         "Ets2"          "Pygb"          "Fgf12"        
    ##  [17] "Me2"           "Tfrc"          "Adcy3"         "Ramp3"        
    ##  [21] "Cops7a"        "Snrpb"         "Retsat"        "Reep1"        
    ##  [25] "Nap1l5"        "Rufy3"         "Psmg1"         "Dclk1"        
    ##  [29] "Carf"          "Vsnl1"         "Sec23a"        "Liph"         
    ##  [33] "Pias2"         "Kif9"          "Hpcal1"        "Med14"        
    ##  [37] "Lmo4"          "Ank"           "Col12a1"       "Lrrc7"        
    ##  [41] "Adam23"        "Lzic"          "Fgf14"         "Tollip"       
    ##  [45] "Gnal"          "Lpin2"         "Rsrc1"         "Necap1"       
    ##  [49] "Atp6v1b2"      "Kdm5c"         "Coq6"          "Rhbdd2"       
    ##  [53] "Mef2c"         "Gabrd"         "Ajap1"         "Pcbp2"        
    ##  [57] "Brwd1"         "Atp6v1g1"      "Mysm1"         "H2afy"        
    ##  [61] "Parp2"         "Fam175b"       "Mad2l1bp"      "Aak1"         
    ##  [65] "Uqcrc2"        "Rragb"         "Cd200"         "Phtf1"        
    ##  [69] "Psd3"          "Rer1"          "Grin2a"        "Atp6ap2"      
    ##  [73] "Slc27a2"       "Ddx19a"        "Cirbp"         "Camk1g"       
    ##  [77] "Srp72"         "Zfp697"        "Cdc5l"         "Prkag2"       
    ##  [81] "Pdk2"          "Gng3"          "Ddx10"         "Farsa"        
    ##  [85] "Nrxn3"         "Hnrnpd"        "Pgs1"          "Fgfr1op"      
    ##  [89] "Nrn1"          "Commd4"        "Msi2"          "Mapk10"       
    ##  [93] "Npc1"          "Cap1"          "Pdss1"         "Eif4b"        
    ##  [97] "Pgk1"          "Stx18"         "Slc20a1"       "Gprasp1"      
    ## [101] "Atp5b"         "Ppp2ca"        "Zcchc17"       "Tubgcp4"      
    ## [105] "Dlat"          "Ypel1"         "Rab2a"         "Gns"          
    ## [109] "Snip1"         "Thumpd3"       "Pcyox1l"       "Pgm2l1"       
    ## [113] "Gdap1"         "Nat10"         "Arih2"         "Prkab1"       
    ## [117] "Rgs4"          "Asap1"         "Amph"          "Prep"         
    ## [121] "Kcnq3"         "Sema3a"        "G6pc3"         "Wdr1"         
    ## [125] "Nek3"          "Ndufa10"       "Vamp1"         "Prkcd"        
    ## [129] "Fosb"          "Mtpap"         "Lgals8"        "Rnf41"        
    ## [133] "Igf1"          "Krba1"         "Pfkm"          "Slc35a3"      
    ## [137] "Stmn2"         "Scg3"          "Psma3"         "Snx17"        
    ## [141] "Ephx2"         "Atp8a2"        "Rfc2"          "Crls1"        
    ## [145] "Daglb"         "Gabra4"        "Dennd3"        "Jak3"         
    ## [149] "Ogfod1"        "Pex1"          "Tmem59"        "Rbm3"         
    ## [153] "Sgsh"          "Bcl2l13"       "Sdc3"          "Stat4"        
    ## [157] "Mlx"           "Kctd5"         "Fads3"         "Trappc2l"     
    ## [161] "Tbc1d9"        "Ptk2b"         "Thnsl2"        "Txn2"         
    ## [165] "Mcph1"         "Snx10"         "Hs3st2"        "Camk2g"       
    ## [169] "Npat"          "Kdelr2"        "Xpnpep1"       "Axin2"        
    ## [173] "Mak16"         "Cacnb2"        "Prpf4"         "Zfp398"       
    ## [177] "Wars"          "Cep57"         "Eif5"          "Ndfip1"       
    ## [181] "Ndrg4"         "Nol6"          "Letmd1"        "Prss23"       
    ## [185] "Gaa"           "Bicd2"         "Dhrs7b"        "M6pr"         
    ## [189] "Usp11"         "Cdc37l1"       "Map2k5"        "Ywhaz"        
    ## [193] "Fen1"          "Atp6ap1"       "Nedd4l"        "Mcam"         
    ## [197] "Vps53"         "Slc17a7"       "Cobll1"        "Actr1a"       
    ## [201] "Tiam2"         "Mrps18b"       "Slc38a7"       "Gabrg2"       
    ## [205] "Pex3"          "Pdcd4"         "4930430F08Rik" "Lrrtm2"       
    ## [209] "Cadps2"        "Mybbp1a"       "Ttc19"         "Lrrfip1"      
    ## [213] "Susd4"         "Mat2b"         "Cyb561"        "Snap47"       
    ## [217] "Sssca1"        "Tomm20"       
    ## [1] "astrocytes_ependymal"
    ## [1] 1
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 1
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 0
    ## [1] "Fold enrichment: 1.17311897842167"
    ## [1] "Standard deviations from mean: 4.51848553378547"
    ## [1] ""
    ## [1] "microglia"
    ## [1] 1
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 1
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 0
    ## [1] "Fold enrichment: 1.13989101452872"
    ## [1] "Standard deviations from mean: 3.68275675202424"
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0
    ## [1] "Fold enrichment: 1.20195269979793"
    ## [1] "Standard deviations from mean: 5.59759745678845"
    ## [1] ""

The results of this analysis can again be plotted using the `ewce.plot` function.

``` r
ewce.plot(tt_results$joint_results)
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-22-1.png)

As was reported in our paper, neuronal genes are found to be downregulated while glial genes are upregulated.

Merging multiple transcriptome studies
--------------------------------------

Where multiple transcriptomic studies have been performed with the same purpose, i.e. seeking differential expression in dlPFC of post-mortem schizophrenics, it is common to want to determine whether they exhibit any shared signal. EWCE can be used to merge the results of multiple studies.

To demonstrate this we use a two further Alzheimer's transcriptome dataset coming from Brodmann areas 36 and 44: these area stored in `tt_alzh_BA36` and `tt_alzh_BA44`. The first step is to run EWCE on each of these individually and store the output into one list.

``` r
# Load the data
data(tt_alzh_BA36)
data(tt_alzh_BA44)

# Run EWCE analysis
tt_results_36 = ewce_expression_data(sct_data=ctd,tt=tt_alzh_BA36,annotLevel=1,ttSpecies="human",sctSpecies="mouse")
```

    ##   [1] "Akap14"        "Arsb"          "Pigv"          "Zzz3"         
    ##   [5] "Jtb"           "Slc39a1"       "Gpm6b"         "Heph"         
    ##   [9] "Rps11"         "Cobl"          "Dlg5"          "Mmaa"         
    ##  [13] "Vgll4"         "Pno1"          "Aga"           "Llgl1"        
    ##  [17] "Sh3glb1"       "Arnt"          "Tfrc"          "Gna11"        
    ##  [21] "Dpf3"          "Thbd"          "Rbm45"         "Rbm15b"       
    ##  [25] "Acadvl"        "Fam107a"       "Ablim2"        "Brd7"         
    ##  [29] "Ccdc77"        "Bard1"         "Smc6"          "Ripk1"        
    ##  [33] "Gzf1"          "Slc35d1"       "Gprc5b"        "Dsg2"         
    ##  [37] "Srprb"         "Prr14"         "Pibf1"         "Znrf2"        
    ##  [41] "Hps3"          "Gabbr1"        "Slc6a19"       "Hkdc1"        
    ##  [45] "Itih5"         "Ikzf4"         "Nkd1"          "Kcnn2"        
    ##  [49] "Cdca7"         "Gpc4"          "Homer3"        "Fam149b"      
    ##  [53] "Zbtb43"        "Slc7a11"       "Smchd1"        "Col27a1"      
    ##  [57] "Elf2"          "Kctd7"         "Taok2"         "Setd3"        
    ##  [61] "Gne"           "Ercc5"         "Efs"           "Ddit4l"       
    ##  [65] "Rap2c"         "Sorbs1"        "Ncl"           "Klhl14"       
    ##  [69] "Snta1"         "Vps37a"        "Gatm"          "Zfp36l2"      
    ##  [73] "Rer1"          "Ubiad1"        "Htra3"         "Htra2"        
    ##  [77] "Prickle1"      "Antxr2"        "Fcrls"         "Lgals1"       
    ##  [81] "Ddx17"         "Masp1"         "Slc38a1"       "Lrig1"        
    ##  [85] "Ccdc97"        "Nhsl1"         "Zc3h10"        "Plcd3"        
    ##  [89] "Rabep1"        "Bag1"          "Slc48a1"       "Rab26"        
    ##  [93] "Derl2"         "Sema5a"        "Mbd6"          "Itch"         
    ##  [97] "Rrad"          "5730480H06Rik" "Ino80d"        "Wdr91"        
    ## [101] "C330018D20Rik" "Acsbg1"        "Add3"          "Myo5c"        
    ## [105] "Rapgef3"       "Kbtbd8"        "Wdr6"          "Gfap"         
    ## [109] "Clu"           "Slc9a3r1"      "Snx5"          "Pik3c3"       
    ## [113] "Plod1"         "Ddx52"         "Clk4"          "Art3"         
    ## [117] "Diaph1"        "Fstl1"         "Litaf"         "Hk2"          
    ## [121] "Plxnb2"        "Srd5a3"        "Apaf1"         "Grik1"        
    ## [125] "Pttg1ip"       "Depdc1b"       "Osbpl2"        "Xpa"          
    ## [129] "Tmem62"        "Sfrp1"         "Cthrc1"        "Tor1aip2"     
    ## [133] "Ctsh"          "Cnot2"         "Zhx3"          "Smg5"         
    ## [137] "Tbc1d16"       "Atr"           "Aplnr"         "Ubxn1"        
    ## [141] "Snx1"          "Tpcn2"         "Irak2"         "Stx6"         
    ## [145] "Cyth1"         "Rad18"         "Psph"          "Ahsa2"        
    ## [149] "Bcl2l13"       "Notch2"        "Rcsd1"         "Mthfsd"       
    ## [153] "Arhgap22"      "Trim35"        "Maml2"         "Vcan"         
    ## [157] "Strap"         "Ppm1f"         "Pus10"         "Man1a2"       
    ## [161] "Kdelr2"        "Thrap3"        "Tbc1d15"       "Nacc2"        
    ## [165] "Sigirr"        "Fbxo4"         "Arhgef7"       "Ahnak"        
    ## [169] "Cpne6"         "Arap1"         "Bgn"           "Edrf1"        
    ## [173] "Ikbkb"         "Atp6v0a2"      "Pbxip1"        "Ccdc138"      
    ## [177] "Tmem5"         "Atp1a2"        "0610009O20Rik" "5031439G07Rik"
    ## [181] "Cnot1"         "Trp53bp2"      "Nfatc2"        "Lrp1"         
    ## [185] "Tmf1"          "Poldip3"       "Adat2"         "Sntb1"        
    ## [189] "Ralgps2"       "Prlr"          "Ehd1"          "Ncor1"        
    ## [193] "Itga7"         "Xylt2"         "Prcp"          "Tbc1d9b"      
    ## [197] "Agt"           "Cog2"          "Slc38a11"      "Ncstn"        
    ## [201] "Mtmr3"         "Stxbp4"        "Lims2"         "Ttc7"         
    ## [205] "C330027C09Rik" "Mfsd8"         "Zfp148"        "Kdm6a"        
    ## [209] "Irf5"          "Tmem121"      
    ## [1] "astrocytes_ependymal"
    ## [1] 0
    ## [1] "Fold enrichment: 1.52076061278877"
    ## [1] "Standard deviations from mean: 6.25895876659442"
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0.03
    ## [1] "Fold enrichment: 1.21453998668281"
    ## [1] "Standard deviations from mean: 2.48025122049512"
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 1
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0.01
    ## [1] "Fold enrichment: 1.25898021010376"
    ## [1] "Standard deviations from mean: 2.67727666978668"
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.02
    ## [1] "Fold enrichment: 1.14886776593126"
    ## [1] "Standard deviations from mean: 2.26068462781992"
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 1
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 1
    ## [1] ""
    ##   [1] "Nfatc1"        "Tbc1d5"        "Zfp142"        "Crim1"        
    ##   [5] "Nsun4"         "Atxn3"         "Ankrd17"       "Ntrk3"        
    ##   [9] "Lmln"          "Usp33"         "Zfp641"        "Vgll4"        
    ##  [13] "Dnajb12"       "Klf13"         "Ppm1k"         "Eml4"         
    ##  [17] "Mknk1"         "Eps15l1"       "Secisbp2l"     "Trpv2"        
    ##  [21] "Wbp4"          "Fam13a"        "Egfr"          "Gria2"        
    ##  [25] "Rufy3"         "Klhl22"        "Tmem144"       "Rbbp6"        
    ##  [29] "Frmd3"         "Sox4"          "Rbmx2"         "Col4a3bp"     
    ##  [33] "Fam155a"       "Spsb1"         "Pdpk1"         "Pou2f1"       
    ##  [37] "Rgs12"         "1700030K09Rik" "Dhx36"         "Fzd8"         
    ##  [41] "Slco3a1"       "Rapgef2"       "Fbxl18"        "Coq9"         
    ##  [45] "Ireb2"         "Ergic2"        "Brd2"          "Gtpbp8"       
    ##  [49] "Rnf115"        "Rap2a"         "Zyg11b"        "Sesn2"        
    ##  [53] "Insr"          "Grid2"         "Akap6"         "Edn1"         
    ##  [57] "Zcrb1"         "Tcf25"         "Tmem53"        "Rspry1"       
    ##  [61] "Neurog2"       "Stk4"          "Vcp"           "Erc1"         
    ##  [65] "Trio"          "Smc5"          "Rc3h2"         "Thumpd1"      
    ##  [69] "Rbm25"         "Klf9"          "Nadk"          "Slitrk2"      
    ##  [73] "Znrd1"         "Emilin2"       "Sh3rf1"        "Tox4"         
    ##  [77] "Shmt1"         "Psd3"          "Pold1"         "Armc9"        
    ##  [81] "Snca"          "Klf7"          "Pcdh9"         "Mapt"         
    ##  [85] "Plekha1"       "Ap3d1"         "Ccdc59"        "Cdk9"         
    ##  [89] "Nrxn3"         "Trnau1ap"      "Zfp652"        "Bbx"          
    ##  [93] "Cap1"          "Emx1"          "Anks1b"        "Kcnk3"        
    ##  [97] "Slc7a1"        "Slc20a1"       "Zfp397"        "Adrb3"        
    ## [101] "Larp4"         "Brd4"          "Tmem97"        "Ifnar2"       
    ## [105] "Ttc3"          "Ldb3"          "Cenpl"         "Pik3r1"       
    ## [109] "Smarca4"       "Oprk1"         "Sh2b2"         "Nav1"         
    ## [113] "Xdh"           "Ak2"           "Zfp800"        "Arih2"        
    ## [117] "Fnbp4"         "Fcrl1"         "Fam184b"       "Ppig"         
    ## [121] "Foxk2"         "Ank3"          "Cog7"          "Kif5c"        
    ## [125] "Sema6a"        "Rims2"         "Tcf4"          "Gpr155"       
    ## [129] "Nktr"          "Wiz"           "Mbd4"          "Zfp90"        
    ## [133] "Kcnb1"         "Mtrf1l"        "Mycbp2"        "Gnao1"        
    ## [137] "Tpcn2"         "Frmpd4"        "Caml"          "Nrxn1"        
    ## [141] "Sobp"          "Egr4"          "Srrm2"         "Irs2"         
    ## [145] "Glg1"          "Maf"           "Ncam1"         "Sbf1"         
    ## [149] "Med28"         "Anapc5"        "Ubn2"          "Hic2"         
    ## [153] "Mcph1"         "Magi1"         "Abcc3"         "Srpk2"        
    ## [157] "Epha4"         "Cda"           "Ash1l"         "Nacc2"        
    ## [161] "Prpf40a"       "Ubn1"          "Chd7"          "Taf1d"        
    ## [165] "Fbxw7"         "Eif5"          "Cep68"         "Lyar"         
    ## [169] "Atf6"          "Mecp2"         "Sqle"          "Zfp207"       
    ## [173] "Zfp609"        "Ecel1"         "Fbxo21"        "Kpnb1"        
    ## [177] "Tmtc2"         "Ddx24"         "Hsf1"          "Ficd"         
    ## [181] "Plxna2"        "Tulp4"         "Epb41l4a"      "Igf2r"        
    ## [185] "Ankrd11"       "Nucks1"        "Wnt7b"         "Ttc17"        
    ## [189] "Runx1t1"       "Slc25a37"      "Ccdc63"        "Kif1a"        
    ## [193] "Slc1a2"        "Cd47"          "Arid1b"        "Pwp1"         
    ## [197] "Pus3"          "Atp1b1"        "Lrrtm2"        "Rp9"          
    ## [201] "Ddx6"          "Adam19"        "Mlec"          "Fnbp1"        
    ## [205] "Mtmr3"         "Gas2l3"        "Rsf1"          "Kcnj11"       
    ## [209] "Gnas"          "Tnpo3"         "Snx19"         "Nrp1"         
    ## [1] "astrocytes_ependymal"
    ## [1] 1
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 1
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 0.13
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0.54
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.72
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 0
    ## [1] "Fold enrichment: 1.10510129143161"
    ## [1] "Standard deviations from mean: 2.91770004092504"
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0.02
    ## [1] "Fold enrichment: 1.07021154488851"
    ## [1] "Standard deviations from mean: 2.10269623288181"
    ## [1] ""

``` r
tt_results_44 = ewce_expression_data(sct_data=ctd,tt=tt_alzh_BA44,annotLevel=1,ttSpecies="human",sctSpecies="mouse")
```

    ##   [1] "Creb3"         "Tmem87a"       "Tyk2"          "Wdr31"        
    ##   [5] "Pex5"          "Vwa1"          "Acss1"         "Prmt2"        
    ##   [9] "Htatip2"       "Rgs10"         "Dock1"         "Dennd5a"      
    ##  [13] "Dhx16"         "Mcf2l"         "Flii"          "Cd274"        
    ##  [17] "Cobl"          "Casc3"         "Arpp19"        "Ppfibp2"      
    ##  [21] "Mfsd10"        "Stk36"         "Gtpbp6"        "Pkdcc"        
    ##  [25] "Gab1"          "Kcnip3"        "Rbbp9"         "Wac"          
    ##  [29] "Ranbp3"        "Brd3"          "Klf16"         "Notch1"       
    ##  [33] "Pld1"          "Zfp445"        "Plxnd1"        "Zfpm1"        
    ##  [37] "Aspa"          "Mid1ip1"       "Rps24"         "Bcan"         
    ##  [41] "Arhgef10l"     "Gzf1"          "Mov10"         "Ralbp1"       
    ##  [45] "Fzd8"          "Zfp689"        "Vars2"         "Acin1"        
    ##  [49] "Rnf220"        "Ptp4a2"        "Insr"          "Pard6a"       
    ##  [53] "Homer3"        "Pdlim2"        "Pld4"          "Ddx49"        
    ##  [57] "Ctcf"          "Wdr53"         "Taok2"         "Tekt1"        
    ##  [61] "Sap18"         "Wwox"          "Patz1"         "Ddx23"        
    ##  [65] "Stac3"         "Zfp36l2"       "Pvr"           "Dhx38"        
    ##  [69] "Cstb"          "Tns3"          "Cbr4"          "Map4k4"       
    ##  [73] "Ppp4r1"        "Limch1"        "Arhgef2"       "Stard7"       
    ##  [77] "Wnk1"          "Fam171a1"      "Cd82"          "Mchr1"        
    ##  [81] "Spock3"        "Fcrls"         "Daxx"          "Tcof1"        
    ##  [85] "Hnrnpab"       "Apbb2"         "Lpar1"         "Arrdc2"       
    ##  [89] "Npc1"          "Tspan17"       "Jade2"         "Baz1b"        
    ##  [93] "Bag1"          "Kat5"          "Heca"          "Stk40"        
    ##  [97] "Cxcr4"         "Fbxo18"        "Aspdh"         "Zfp637"       
    ## [101] "Larp6"         "Thap11"        "Mylip"         "Smarcc1"      
    ## [105] "Dusp10"        "Dlg1"          "Ezr"           "Rapgef3"      
    ## [109] "Hexim1"        "Fbl"           "Gigyf2"        "Kbtbd4"       
    ## [113] "Pomt2"         "Pum2"          "Cited4"        "Rnf207"       
    ## [117] "Fuk"           "Ube4b"         "Sap30bp"       "Ankzf1"       
    ## [121] "Gjb1"          "3110079O15Rik" "Git1"          "Ggt7"         
    ## [125] "Fbxw8"         "Tbc1d16"       "Ndufv3"        "Hes6"         
    ## [129] "Gsg1l"         "Ascl1"         "Ric8b"         "Osbpl7"       
    ## [133] "Lrrc8d"        "Nf2"           "Ubac1"         "Fbxw4"        
    ## [137] "Myo10"         "Epb41l2"       "Ahsa2"         "Nfasc"        
    ## [141] "Pusl1"         "Maf"           "Trabd"         "Mthfsd"       
    ## [145] "Smarcd3"       "Ttll9"         "Arhgap22"      "Zfp24"        
    ## [149] "Vcan"          "Safb2"         "Gpr108"        "Crtap"        
    ## [153] "Qdpr"          "Acy3"          "Zfp174"        "Oraov1"       
    ## [157] "Plekhb1"       "Spop"          "Nlgn3"         "Fam57b"       
    ## [161] "Pkp4"          "Tmem184b"      "Exoc3"         "Shroom1"      
    ## [165] "Tmcc3"         "Trim41"        "Ankrd13b"      "Fam173b"      
    ## [169] "Mex3a"         "Rasgrp3"       "Ptpn11"        "Exosc3"       
    ## [173] "Sqle"          "Gtf3c5"        "Slc38a2"       "Hs1bp3"       
    ## [177] "Atp1b3"        "Dusp7"         "Sh3gl3"        "Tmem101"      
    ## [181] "Aig1"          "Tmcc2"         "Nrip2"         "Bad"          
    ## [185] "Nfkbib"        "Tmem119"       "Adamts13"      "Rnf130"       
    ## [189] "Ppm1d"         "Cdr2l"         "Ncor1"         "Srrt"         
    ## [193] "Mcm7"          "Vsig2"         "Cd84"          "Sccpdh"       
    ## [197] "Selplg"        "Capn5"         "Itpk1"         "Rcor2"        
    ## [201] "Pgpep1"        "R3hcc1"        "Chka"          "Vezf1"        
    ## [205] "Nap1l1"        "Fbxo32"       
    ## [1] "astrocytes_ependymal"
    ## [1] 0.72
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0.79
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 1
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0.11
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0
    ## [1] "Fold enrichment: 1.48069228148824"
    ## [1] "Standard deviations from mean: 7.77128885133655"
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 1
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0.92
    ## [1] ""
    ##   [1] "Lgals3"        "Rab36"         "Uba2"          "Mab21l2"      
    ##   [5] "Dock7"         "Vegfc"         "Fez2"          "Stk17b"       
    ##   [9] "Ica1"          "Bckdhb"        "Ntrk3"         "Nqo1"         
    ##  [13] "Aim2"          "Dusp11"        "Ttc21b"        "Ahcyl1"       
    ##  [17] "Arrdc3"        "Pno1"          "Zscan12"       "Sardh"        
    ##  [21] "Ppm1k"         "Rpl35a"        "Eps15l1"       "Cdh8"         
    ##  [25] "Stk38l"        "Pla2r1"        "Rragc"         "Tph1"         
    ##  [29] "Cmpk2"         "Wwp2"          "Cpne8"         "Cdca3"        
    ##  [33] "Slamf8"        "Det1"          "Mki67"         "Med31"        
    ##  [37] "Fgf7"          "Tubb2b"        "Nup214"        "Cd38"         
    ##  [41] "Fosl2"         "Mlf1"          "Ola1"          "Pigo"         
    ##  [45] "Rarres1"       "Ttr"           "Wrap53"        "Srd5a1"       
    ##  [49] "Gbp5"          "Ccdc68"        "Sv2c"          "Pdk1"         
    ##  [53] "Adcyap1"       "Chmp1b"        "Gadd45a"       "Bicc1"        
    ##  [57] "Otud4"         "Luzp2"         "Mrps30"        "Gls"          
    ##  [61] "Pou3f3"        "Itga8"         "Dync2li1"      "Rlim"         
    ##  [65] "Mmp24"         "Cd40"          "Fkbp14"        "Ide"          
    ##  [69] "Setdb2"        "Abhd2"         "Il6st"         "Dnase1l3"     
    ##  [73] "Cds1"          "Ptprg"         "Cp"            "Tuft1"        
    ##  [77] "Lmcd1"         "Gabrr1"        "Rfesd"         "Slc35e1"      
    ##  [81] "Dsc3"          "Ddx58"         "Cyb5a"         "Hectd2"       
    ##  [85] "Psd3"          "Cstf3"         "Sec22a"        "Tmtc1"        
    ##  [89] "Pdia3"         "Dph5"          "Car12"         "Atl2"         
    ##  [93] "Efemp1"        "Mchr1"         "Atp5g2"        "Smad3"        
    ##  [97] "Hnrnpd"        "Slc22a4"       "Aste1"         "Slc16a6"      
    ## [101] "Cenpj"         "Tmem97"        "Phtf2"         "Fmnl3"        
    ## [105] "Nat10"         "Hsf4"          "Agfg1"         "Rnh1"         
    ## [109] "Mllt3"         "Bcap29"        "Dnajc18"       "Prkab1"       
    ## [113] "Hmgcl"         "Nid1"          "Stx17"         "Fbxw2"        
    ## [117] "Nudcd1"        "Hddc2"         "Kcnq3"         "Smyd4"        
    ## [121] "Ppp1r3d"       "Spcs2"         "Lgals8"        "Ddx31"        
    ## [125] "Ppp2r5e"       "Dtna"          "Sf1"           "Tnfaip2"      
    ## [129] "Rfc5"          "8430408G22Rik" "Tlr1"          "Sdc1"         
    ## [133] "Ccdc96"        "Ogfod1"        "Pex1"          "Rbm3"         
    ## [137] "Dusp15"        "Ahsa2"         "Sgsh"          "Arhgap24"     
    ## [141] "Mboat7"        "Ttc26"         "Pacrg"         "AU019823"     
    ## [145] "Nip7"          "Lsm11"         "Bmpr2"         "Ctnnal1"      
    ## [149] "Zw10"          "Slc31a1"       "St7l"          "Ezh1"         
    ## [153] "Zfr2"          "Adam22"        "Mpzl1"         "Etnk1"        
    ## [157] "Rpl27a"        "Itsn1"         "Zbtb25"        "Apool"        
    ## [161] "Gtf3c4"        "Prss23"        "Elavl1"        "Pop1"         
    ## [165] "Sgsm1"         "Ublcp1"        "Usp15"         "Chpt1"        
    ## [169] "Fbln1"         "Slit3"         "Kcnip1"        "Hgf"          
    ## [173] "Cd84"          "Sft2d2"        "Slc30a10"      "Rasal2"       
    ## [177] "Susd4"         "Rsf1"          "Sparcl1"       "Gabrb2"       
    ## [181] "Rbck1"        
    ## [1] "astrocytes_ependymal"
    ## [1] 0.04
    ## [1] "Fold enrichment: 1.16434783921907"
    ## [1] "Standard deviations from mean: 1.68751807472873"
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0.09
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 0.97
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0
    ## [1] "Fold enrichment: 1.25546279672692"
    ## [1] "Standard deviations from mean: 2.71084273786502"
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.82
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 0.89
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0.98
    ## [1] ""

``` r
# Fill a list with the results
results = add.res.to.merging.list(tt_results)
results = add.res.to.merging.list(tt_results_36,results)
results = add.res.to.merging.list(tt_results_44,results)

# Perform the merged analysis
merged_res = merged_ewce(results,reps=10) # <- For publication reps should be higher
print(merged_res)
```

    ##                                   CellType       p        fc sd_from_mean
    ## astrocytes_ependymal  astrocytes_ependymal 0.00000 1.2890874    6.0667151
    ## endothelial-mural        endothelial-mural 0.00019 1.1767707    3.8086661
    ## interneurons                  interneurons 1.00000 0.8350048   -7.2614211
    ## microglia                        microglia 0.00001 1.2412430    4.1893293
    ## oligodendrocytes          oligodendrocytes 0.00000 1.2229410    6.1150611
    ## pyramidal CA1                pyramidal CA1 1.00000 0.8632602   -6.3231402
    ## pyramidal SS                  pyramidal SS 1.00000 0.8771292   -5.8447364
    ## astrocytes_ependymal1 astrocytes_ependymal 0.99659 0.8751966   -2.6584406
    ## endothelial-mural1       endothelial-mural 0.99903 0.8642894   -2.7343326
    ## interneurons1                 interneurons 0.01058 1.0538325    2.2561299
    ## microglia1                       microglia 0.81375 0.9450484   -0.8981942
    ## oligodendrocytes1         oligodendrocytes 0.99974 0.8805943   -2.8813514
    ## pyramidal CA11               pyramidal CA1 0.00026 1.0719668    3.2597619
    ## pyramidal SS1                 pyramidal SS 0.00023 1.0762653    3.6619888
    ##                       Direction
    ## astrocytes_ependymal         Up
    ## endothelial-mural            Up
    ## interneurons                 Up
    ## microglia                    Up
    ## oligodendrocytes             Up
    ## pyramidal CA1                Up
    ## pyramidal SS                 Up
    ## astrocytes_ependymal1      Down
    ## endothelial-mural1         Down
    ## interneurons1              Down
    ## microglia1                 Down
    ## oligodendrocytes1          Down
    ## pyramidal CA11             Down
    ## pyramidal SS1              Down

The results can then be plotted as normal using the `ewce.plot` function.

``` r
print(ewce.plot(merged_res))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-24-1.png)

The merged results from all three Alzheimer's brain regions are found to be remarkably similar, as was reported in our paper.

Conditional cell type enrichment analysis
=========================================

Controlling for expression in another cell type
-----------------------------------------------

In a followup paper we found that an enrichment detected for Schizophrenia in Somatosensory Pyramidal neurons could be explained by accounting for expression in Hippocampal CA1 pyramidal neurons. These results are described here:

    Nathan G. Skene, Julien Bryois, Trygve E. Bakken, Gerome Breen, James J. Crowley, Helena Gaspar, Paola Giusti-Rodriguez, Rebecca D. Hodge, Jeremy A. Miller, Ana Munoz-Manchado, Michael C. O'Donovan, Michael J. Owen, Antonio F. Pardinas, Jesper Ryge, James T. R. Walters, Sten Linnarsson, Ed S. Lein, - Major Depressive Disorder Working Group of the PGC, Patrick F. Sullivan, Jens Hjerling-Leffler
    "Genetic Identification Of Brain Cell Types Underlying Schizophrenia."
    bioRxiv, 2017

Those results were generated using an alternative enrichment method designed for use with GWAS Summary Statistics rather than gene sets. The same sort of approach can be extended to EWCE as well, and we have implemented it within this package. When testing for enrichment the other gene sets that are sampled are selected to have equivilent specificity in the controlled celltype.

We demonstrate it's use below to test whether the enrichment in astrocytes is still present after controlling for the enrichment within microglia:

``` r
data("mouse_to_human_homologs")
m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
mouse.hits = unique(m2h[m2h$HGNC.symbol %in% example_genelist,"MGI.symbol"])
mouse.bg  = unique(m2h$MGI.symbol)

reps=1000
unconditional_results = bootstrap.enrichment.test(sct_data=ctd,hits=mouse.hits,
                    bg=mouse.bg,reps=reps,annotLevel=1,genelistSpecies="mouse",sctSpecies="mouse")
```

    ##  [1] "Apoe"    "Inpp5d"  "Cd2ap"   "Nme8"    "Cass4"   "Mef2c"   "Zcwpw1" 
    ##  [8] "Bin1"    "Clu"     "Celf1"   "Abca7"   "Slc24a4" "Ptk2b"   "Picalm" 
    ## [15] "Fermt2"  "Sorl1"  
    ## [1] "astrocytes_ependymal"
    ## [1] 0.044
    ## [1] "Fold enrichment: 1.61075434352815"
    ## [1] "Standard deviations from mean: 1.90933629115277"
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0.604
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 1
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0.015
    ## [1] "Fold enrichment: 2.06505599892457"
    ## [1] "Standard deviations from mean: 2.7209383544878"
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.522
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 0.658
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0.695
    ## [1] ""

``` r
conditional_results_micro = bootstrap.enrichment.test(sct_data=ctd,hits=mouse.hits,
                    bg=mouse.bg,reps=reps,annotLevel=1,controlledCT="microglia",genelistSpecies="mouse",sctSpecies="mouse")
```

    ##  [1] "Apoe"    "Inpp5d"  "Cd2ap"   "Nme8"    "Cass4"   "Mef2c"   "Zcwpw1" 
    ##  [8] "Bin1"    "Clu"     "Celf1"   "Abca7"   "Slc24a4" "Ptk2b"   "Picalm" 
    ## [15] "Fermt2"  "Sorl1"  
    ## [1] "astrocytes_ependymal"
    ## [1] 0.022
    ## [1] "Fold enrichment: 1.71342601619128"
    ## [1] "Standard deviations from mean: 2.44919592202938"
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0.421
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 1
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0.535
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.399
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 0.4
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0.457
    ## [1] ""

``` r
conditional_results_astro = bootstrap.enrichment.test(sct_data=ctd,hits=mouse.hits,
                    bg=mouse.bg,reps=reps,annotLevel=1,controlledCT="astrocytes_ependymal",genelistSpecies="mouse",sctSpecies="mouse")
```

    ##  [1] "Apoe"    "Inpp5d"  "Cd2ap"   "Nme8"    "Cass4"   "Mef2c"   "Zcwpw1" 
    ##  [8] "Bin1"    "Clu"     "Celf1"   "Abca7"   "Slc24a4" "Ptk2b"   "Picalm" 
    ## [15] "Fermt2"  "Sorl1"  
    ## [1] "astrocytes_ependymal"
    ## [1] 0.193
    ## [1] ""
    ## [1] "endothelial-mural"
    ## [1] 0.389
    ## [1] ""
    ## [1] "interneurons"
    ## [1] 1
    ## [1] ""
    ## [1] "microglia"
    ## [1] 0.005
    ## [1] "Fold enrichment: 2.10197457213456"
    ## [1] "Standard deviations from mean: 2.95302989520271"
    ## [1] ""
    ## [1] "oligodendrocytes"
    ## [1] 0.483
    ## [1] ""
    ## [1] "pyramidal CA1"
    ## [1] 0.487
    ## [1] ""
    ## [1] "pyramidal SS"
    ## [1] 0.524
    ## [1] ""

``` r
full_res1 = data.frame(unconditional_results$results,list="Unconditional Enrichment")
full_res2 = data.frame(conditional_results_micro$results,list="Conditional Enrichment (Microglia controlled)")
full_res3 = data.frame(conditional_results_astro$results,list="Conditional Enrichment (Astrocyte controlled)")
merged_results = rbind(rbind(full_res1,full_res2),full_res3)
print(ewce.plot(total_res=merged_results,mtc_method="BH"))
```

![](/private/var/folders/qk/f3p79d9523j3kjt3qkrcn0k0000csj/T/RtmpXwGs21/preview-ca6770182eef.dir/EWCE_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-25-1.png)

When controlling for astrocytes the enrichment is astrocytes is totally abolished as expected, and vica versa. The enrichment in microglia remains strongly significant however after controlling for microglia, suggesting that this enrichment is independent of that in astrocytes.

Gene set enrichment analysis controlling for cell type expression
-----------------------------------------------------------------

Traditionally the standard analysis run on all gene sets was the GO enrichment analysis. Once you have established that a given gene list is enriched for a given celltype, it becomes questionable whether a GO enrichment is driven purely by the underlying cell type enrichment. For instance, it is well established that genes associated with schizophrenia are enriched for the human post-synaptic density genes, however, it has also been shown that schizophrenia is enriched for specificity in CA1 pyramidal neurons (which highly express hPSD genes). These two enrichments can be disassociated using the following analysis:

``` r
data("mouse_to_human_homologs")
m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])

data("schiz_genes")
data("id_genes")
mouse.hits.schiz = unique(m2h[m2h$HGNC.symbol %in% schiz_genes,"MGI.symbol"])
mouse.hits.id = unique(m2h[m2h$HGNC.symbol %in% id_genes,"MGI.symbol"])
mouse.bg  = unique(m2h$MGI.symbol)

data("hpsd_genes")
mouse.hpsd = unique(m2h[m2h$HGNC.symbol %in% hpsd_genes,"MGI.symbol"])

data("rbfox_genes")

res_hpsd_schiz = controlled_geneset_enrichment(disease_genes=mouse.hits.schiz, functional_genes = mouse.hpsd, bg_genes = mouse.bg, sct_data = ctd, annotLevel = 1, reps=1000, controlledCT="pyramidal CA1")
res_rbfox_schiz = controlled_geneset_enrichment(disease_genes=mouse.hits.schiz, functional_genes = rbfox_genes, bg_genes = mouse.bg, sct_data = ctd, annotLevel = 1, reps=1000, controlledCT="pyramidal CA1")
print(res_hpsd_schiz)
```

    ## $p_controlled
    ## [1] 0.124
    ## 
    ## $z_controlled
    ## [1] 1.049584
    ## 
    ## $p_uncontrolled
    ## [1] 0.046
    ## 
    ## $z_uncontrolled
    ## [1] 1.636406
    ## 
    ## $reps
    ## [1] 1000
    ## 
    ## $controlledCT
    ## [1] "pyramidal CA1"
    ## 
    ## $actualOverlap
    ## [1] 54

``` r
print(res_rbfox_schiz)
```

    ## $p_controlled
    ## [1] 0.001
    ## 
    ## $z_controlled
    ## [1] 3.553008
    ## 
    ## $p_uncontrolled
    ## [1] 0
    ## 
    ## $z_uncontrolled
    ## [1] 4.438384
    ## 
    ## $reps
    ## [1] 1000
    ## 
    ## $controlledCT
    ## [1] "pyramidal CA1"
    ## 
    ## $actualOverlap
    ## [1] 163

``` r
res_hpsd_id = controlled_geneset_enrichment(disease_genes=mouse.hits.id, functional_genes = mouse.hpsd, bg_genes = mouse.bg, sct_data = ctd, annotLevel = 1, reps=1000, controlledCT="pyramidal SS")
res_rbfox_id = controlled_geneset_enrichment(disease_genes=mouse.hits.id, functional_genes = rbfox_genes, bg_genes = mouse.bg, sct_data = ctd, annotLevel = 1, reps=1000, controlledCT="pyramidal SS")
print(res_hpsd_id)
```

    ## $p_controlled
    ## [1] 0
    ## 
    ## $z_controlled
    ## [1] 4.440125
    ## 
    ## $p_uncontrolled
    ## [1] 0
    ## 
    ## $z_uncontrolled
    ## [1] 5.010771
    ## 
    ## $reps
    ## [1] 1000
    ## 
    ## $controlledCT
    ## [1] "pyramidal SS"
    ## 
    ## $actualOverlap
    ## [1] 136

``` r
print(res_rbfox_id)
```

    ## $p_controlled
    ## [1] 0
    ## 
    ## $z_controlled
    ## [1] 3.706117
    ## 
    ## $p_uncontrolled
    ## [1] 0
    ## 
    ## $z_uncontrolled
    ## [1] 3.969491
    ## 
    ## $reps
    ## [1] 1000
    ## 
    ## $controlledCT
    ## [1] "pyramidal SS"
    ## 
    ## $actualOverlap
    ## [1] 310

The analysis also tests for enrichment of Rbfox binding genes in the schizophrenia susceptibility genes, as well as both hPSD and Rbfox genes in Intellectual Disability genes. All of the enrichments are confirmed as still being present after controlling for the associated cell type, apart from the enrichment of PSD genes in Schizophrenia which falls from borderline to non-significant.

Controlling for multiple cell types
-----------------------------------

``` r
controlledCTs = c("pyramidal CA1","pyramidal SS","interneurons")

res_hpsd_schiz = controlled_geneset_enrichment(disease_genes=mouse.hits.schiz, functional_genes = mouse.hpsd, bg_genes = mouse.bg, sct_data = ctd, annotLevel = 1, reps=1000, controlledCT=controlledCTs)
res_rbfox_schiz = controlled_geneset_enrichment(disease_genes=mouse.hits.schiz, functional_genes = rbfox_genes, bg_genes = mouse.bg, sct_data = ctd, annotLevel = 1, reps=1000, controlledCT=controlledCTs)
print(res_hpsd_schiz)
```

    ## $p_controlled
    ## [1] 0.215
    ## 
    ## $z_controlled
    ## [1] 0.7203522
    ## 
    ## $p_uncontrolled
    ## [1] 0.038
    ## 
    ## $z_uncontrolled
    ## [1] 1.665142
    ## 
    ## $reps
    ## [1] 1000
    ## 
    ## $controlledCT
    ## [1] "pyramidal CA1" "pyramidal SS"  "interneurons" 
    ## 
    ## $actualOverlap
    ## [1] 54

``` r
print(res_rbfox_schiz)
```

    ## $p_controlled
    ## [1] 0
    ## 
    ## $z_controlled
    ## [1] 3.41526
    ## 
    ## $p_uncontrolled
    ## [1] 0
    ## 
    ## $z_uncontrolled
    ## [1] 4.492349
    ## 
    ## $reps
    ## [1] 1000
    ## 
    ## $controlledCT
    ## [1] "pyramidal CA1" "pyramidal SS"  "interneurons" 
    ## 
    ## $actualOverlap
    ## [1] 163

``` r
res_hpsd_id = controlled_geneset_enrichment(disease_genes=mouse.hits.id, functional_genes = mouse.hpsd, bg_genes = mouse.bg, sct_data = ctd, annotLevel = 1, reps=1000, controlledCT=controlledCTs)
res_rbfox_id = controlled_geneset_enrichment(disease_genes=mouse.hits.id, functional_genes = rbfox_genes, bg_genes = mouse.bg, sct_data = ctd, annotLevel = 1, reps=1000, controlledCT=controlledCTs)
print(res_hpsd_id)
```

    ## $p_controlled
    ## [1] 0
    ## 
    ## $z_controlled
    ## [1] 4.930243
    ## 
    ## $p_uncontrolled
    ## [1] 0
    ## 
    ## $z_uncontrolled
    ## [1] 4.839632
    ## 
    ## $reps
    ## [1] 1000
    ## 
    ## $controlledCT
    ## [1] "pyramidal CA1" "pyramidal SS"  "interneurons" 
    ## 
    ## $actualOverlap
    ## [1] 136

``` r
print(res_rbfox_id)
```

    ## $p_controlled
    ## [1] 0
    ## 
    ## $z_controlled
    ## [1] 4.142594
    ## 
    ## $p_uncontrolled
    ## [1] 0
    ## 
    ## $z_uncontrolled
    ## [1] 3.784632
    ## 
    ## $reps
    ## [1] 1000
    ## 
    ## $controlledCT
    ## [1] "pyramidal CA1" "pyramidal SS"  "interneurons" 
    ## 
    ## $actualOverlap
    ## [1] 310

References
==========

1. Skene, N. & Grant, S. Identification of vulnerable cell types in major brain disorders using single cell transcriptomes and expression weighted cell type enrichment. *Frontiers in Neuroscience* (2016). doi:[10.3389/fnins.2016.00016](https://doi.org/10.3389/fnins.2016.00016)

2. Zeisel, A. *et al.* Cell types in the mouse cortex and hippocampus revealed by single-cell rna-seq. *Science* **347,** 1138–1142 (2015).

3. Haroutunian, V., Katsel, P. & Schmeidler, J. Transcriptional vulnerability of brain regions in alzheimer’s disease and dementia. *Neurobiology of aging* **30,** 561–573 (2009).
