# test merge_two_expfiles
test_that("merging two expression files", {
  #use vignette data
  
  data(cortex_mrna)
  #NOTE: Can't use hypothalamus data for test as issue with travis and read_xl::read_excel()
  
  
  # Download the hypothalamus data and unzip
  #if(!file.exists("GSE74672_expressed_mols_with_classes.xlsx")){
  #  download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74672/suppl/GSE74672_expressed_mols_with_classes.xlsx.gz", destfile="GSE74672_expressed_mols_with_classes.xlsx.gz")
  #  system("gunzip GSE74672_expressed_mols_with_classes.xlsx.gz")
  #}
  
  # Read in the hypothalamus data
  #library(readxl)
  #hypo_dat = read_excel("GSE74672_expressed_mols_with_classes.xlsx")
  
  # Extract the expression data, gene symbols and annotation data
  #exp = data.matrix(hypo_dat[12:dim(hypo_dat)[1],2:dim(hypo_dat)[2]])
  #rownames(exp) = data.frame(hypo_dat[12:dim(hypo_dat)[1],1])[,1]
  #level1class = data.frame(level1class=t(hypo_dat[1,2:dim(hypo_dat)[2]]),stringsAsFactors = FALSE)[,1]
  #level2class = data.frame(leve2class=t(hypo_dat[2,2:dim(hypo_dat)[2]]),stringsAsFactors = FALSE)[,1]
  #cell_id     = colnames(hypo_dat)[2:dim(hypo_dat)[2]]
  #hypo_annot  = data.frame(cell_id=cell_id,level1class=level1class,level2class=level2class,stringsAsFactors = FALSE)
  
  # Drop the glia and unclassified cells (which don't have level 2  annotations)
 # hypo_annot  = hypo_annot[!is.na(hypo_annot$level2class) & !hypo_annot$level2class=="uc",]
 # hypo_exp    = exp[,hypo_annot$cell_id]
  
  # Make the celltype names more aesthetically pleasing
 # hypo_annot$level2class=gsub(",",";",hypo_annot$level2class)
#  hypo_annot$level1class[grep("Oxt;|^Avp",hypo_annot$level2class)] = "Oxytocin / Vasopressin Expressing Neurons"
#  hypo_annot$level1class[grep("^Th;|^Dopamine",hypo_annot$level2class)] = "Hypothalamic Dopaminergic Neurons"
#  hypo_annot$level1class[grepl("^Vglut2|^Trh|^Qrfp|^Hcrt|^Pmch|^Adcyap1|^Npvf|^Ghrh|^Hmit|^Nms|^Vip;|^Per2|Tnr$|^Gad-low;Gnrh",hypo_annot$level2class) & grepl("neurons",hypo_annot$level1class)] = "Hypothalamic Glutamatergic Neurons"
#  hypo_annot$level1class[grepl("GABA|^Sst|^Crh|^Npy|^Pomc|^Galanin|^Otof|Pnoc$|^Calcr-high",hypo_annot$level2class) & grepl("^neurons$",hypo_annot$level1class)] = "Hypothalamic GABAergic Neurons"
#  hypo_annot$level2class[hypo_annot$level2class!=""] = sprintf("Hypothalamic %s Neuron",hypo_annot$level2class[hypo_annot$level2class!=""])
  
  # Fix bad MGI symbols
  #hypo_exp_CORRECTED = fix.bad.mgi.symbols(hypo_exp)
  # Merge the datasets
  merged_KI = EWCE::merge_two_expfiles(exp1=cortex_mrna$exp,#hypo_exp,#hypo_exp_CORRECTED,  
                                 exp2=cortex_mrna$exp,
                                 annot1=cortex_mrna$annot,#hypo_annot,        
                                 annot2=cortex_mrna$annot,
                                 name1="Hypothalamus (KI)", name2="Cortex/Hippo (KI)")
  
  #check merge didn't drop any samples
  test1 <- all.equal(sort(colnames(merged_KI$exp)),
                     sort(c(colnames(cortex_mrna$exp),
                            colnames(cortex_mrna$exp))))
  #remove folder once tested
  #unlink("GSE74672_expressed_mols_with_classes.xlsx")
  
  #fail if either function didn't create directory and add files
  expect_equal(test1,TRUE)
})
