#' Table of Human-->Mouse orthologs for all human genes
#'
#' A dataset containing the MGI and HGNC symbols, Human and Mouse Entrez and Ensembl gene IDs for all human
#' orthologs for mouse genes. Whenin the mouse genes are defined based on a list of all MGI markers from the
#' MGI website (downloaded as MRK_List2.rpt file from http://goo.gl/mjf2GQ)
#'
#' @source
#' The code to prepare the .Rda file file from the marker file is:
#' \code{
#' markers = read.csv("MRK_List2.rpt",sep="\t")
#' genes = markers[markers$Feature.Type=="protein coding gene",]
#' listMarts(host="www.ensembl.org")
#' human <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
#' mouse <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
#' mouse_to_human_homologs = getLDS(attributes = c("mgi_symbol","entrezgene","ensembl_gene_id"),
#'                                  filters = "mgi_symbol", values = genes$Marker.Symbol,
#'                                  mart = mouse,
#'                                  attributesL = c("hgnc_symbol","ensembl_gene_id","entrezgene"), martL = human)
#' unique_mgi = setdiff(mouse_to_human_homologs$MGI.symbol,mouse_to_human_homologs$MGI.symbol[duplicated(mouse_to_human_homologs$MGI.symbol)])
#' unique_hgnc = setdiff(mouse_to_human_homologs$HGNC.symbol,mouse_to_human_homologs$HGNC.symbol[duplicated(mouse_to_human_homologs$HGNC.symbol)])
#' mouse_to_human_homologs = mouse_to_human_homologs[mouse_to_human_homologs$MGI.symbol %in% unique_mgi & mouse_to_human_homologs$HGNC.symbol %in% unique_hgnc,]
#' save(mouse_to_human_homologs,file="mouse_to_human_homologs.Rda")
#' }
#'
"mouse_to_human_homologs"

#' All MGI gene symbols
#'
#' A dataset containing all MGI symbols
#'
#' @source
#' The code to prepare the .Rda file file from the marker file is:
#' \code{
#' library("biomaRt")
#' listMarts(host="www.ensembl.org")
#' mouse <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
#' mgi_symbols = getBM(attributes=c("mgi_symbol","ensembl_gene_id"), mart=mouse)
#' all_mgi = unique(mgi_symbols$mgi_symbol)
#' save(all_mgi,file="all_mgi.Rda")
#' }
#'
"all_mgi"

#' All HGNC gene symbols
#'
#' A dataset containing all HGNC symbols
#'
#' @source
#' The code to prepare the .Rda file file from the marker file is:
#' \code{
#' library("biomaRt")
#' listMarts(host="www.ensembl.org")
#' human <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
#' hgnc_symbols = getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), mart=human)
#' all_hgnc = unique(hgnc_symbols$hgnc_symbol)
#' save(all_hgnc,file="all_hgnc.Rda")
#' }
#'
"all_hgnc"

#' Example Gene list
#'
#' A list of genes genetically associated with Alzheimer's disease.
#' @source
#' These were obtained from two sources:
#' http://www.alzgene.org/TopResults.asp and PMID24162737
"example_genelist"

#' Example table of differential expression (Alzheimer's BA46)
#'
#' A list of genes found to be differentially expressed in the BA46 in Alzheimer's disease.
#' @source
#' The table was determined based on data associated with the paper with PMID:17845826
"tt_alzh"

#' Example table of differential expression (Alzheimer's BA36)
#'
#' A list of genes found to be differentially expressed in the BA36 in Alzheimer's disease.
#' @source
#' The table was determined based on data associated with the paper with PMID:17845826
"tt_alzh_BA36"

#' Example table of differential expression (Alzheimer's BA44)
#'
#' A list of genes found to be differentially expressed in the BA44 in Alzheimer's disease.
#' @source
#' The table was determined based on data associated with the paper with PMID:17845826
"tt_alzh_BA44"

#' The genes from Karolinska cortex/hippocampus and hypothalamus single cell transcriptome datasets
#'
#' All genes from the merged SCT dataset
#' @source
#' The datasets were downloaded from the website associated with the paper & GEO and loaded using read_celltype_data PMID:25700174
"ctd"

#' Schizophrenia susceptibility genes from CLOZUK
#' @source
#' Extended data tables downloaded from: http://www.biorxiv.org/content/early/2016/08/09/068593.figures-only
#' schiz_genes = read.csv("/Users/ns9/Google Drive/DiseaseLists/Schizophrenia_CLOZUK_geneWideSignif.txt",stringsAsFactors = FALSE)[-1,1]
#' save(schiz_genes,file="data/schiz_genes.rda")
'schiz_genes'

#' Human post-synaptic density
#' @source
#' From Supplementary Table 2 of https://www.nature.com/neuro/journal/v14/n1/full/nn.2719.html
#' hpsd_genes = read.csv("/Users/ns9/Google Drive/DiseaseLists/hPSD.txt",stringsAsFactors = FALSE)[-1,1]
#' save(hpsd_genes,file="data/hpsd_genes.rda")
'hpsd_genes'

#' Rbfox binding genes
#' @source
#' from supplementary table 1 of HITS-CLIP and Integrative Modeling Define the Rbfox Splicing-Regulatory Network Linked to Brain Development and Autism. All with rbfox2 count greater than 4 or summed rbfox 1 and 3 greater than 12
#' rbfox_genes = read.csv("/Users/ns9/Google Drive/DiseaseLists/Rbfox_binding.txt",stringsAsFactors = FALSE)[-1,1]
#' save(rbfox_genes,file="data/rbfox_genes.rda")
'rbfox_genes'

#' Intellectual disability genes
#' @source
#' from http://compbio.charite.de/hpoweb/showterm?id=HP:0001249
#' id_genes = read.csv("/Users/ns9/Google Drive/DiseaseLists/Intellectual Disability July2017.txt",stringsAsFactors = FALSE)[-1,1]
#' save(id_genes,file="data/id_genes.rda")
'id_genes'

#' Karolinska Cortex/Hippocamus dataset subsample
#' @source
#' download.file("goo.gl/r5Y24y",destfile="expression_mRNA_17-Aug-2014.txt")
#' path = "expression_mRNA_17-Aug-2014.txt"
#' cortex_mrna  = load.linnarsson.sct.data("expression_mRNA_17-Aug-2014.txt")
#' save(cortex_mrna,file="data/cortex_mrna.rda")
'cortex_mrna'

#' Alzheimers disease top100 GWAS genes generated by MAGMA using the iGAP summary statistics
#' @source
#' hgnc_symbols = getBM(attributes=c("hgnc_symbol","ensembl_gene_id","entrezgene"), mart=human)
#' magma2=merge(magma,hgnc_symbols,by="GENE")
#' write.csv(magma2[1:100,]$hgnc_symbol,file="Alzh_IGAP_top100magma.txt",quote=FALSE,row.names = FALSE)
'alzh_gwas_top100'

#' MGI synonym data
#' @source
#' download.file("http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt", destfile="MRK_List2.rpt")
#' mgi_synonym_data = read.csv(mrk_file_path,sep="\t",stringsAsFactors = FALSE)
#' mgi_synonym_data = mgi_data[!mgi_data$Marker.Synonyms..pipe.separated.=="",]
#' save(mgi_synonym_data,file="data/mgi_synonym_data.rda", compress='xz')
'mgi_synonym_data'
