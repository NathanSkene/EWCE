# test load.linnarsson.sct.data
test_that("loading of linnarsson dataset", {
  #from vignette
  
  download.file("goo.gl/r5Y24y",
                destfile="expression_mRNA_17-Aug-2014.txt") 
  
  path = "expression_mRNA_17-Aug-2014.txt"
  
  cortex_mrna  = load.linnarsson.sct.data(path)
  
  #remove folder once tested
  unlink("expression_mRNA_17-Aug-2014.txt")
  
  #fail if data didn't load correctly
  expect_equal(
    all(length(cortex_mrna)==2 &&
        class(cortex_mrna$exp) =="matrix" && 
        nrow(cortex_mrna$exp)>1 &&
        class(cortex_mrna$annot)=="data.frame" && 
        nrow(cortex_mrna$annot)>1),
  TRUE)
  
})