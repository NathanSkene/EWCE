#Test for fix.bad.mgi.symbols
test_that("method to remove/fix an expected set of genes", {
  #Use Vignette Dataset to check function, alter input gene names
  data(cortex_mrna)
  if(!file.exists("MRK_List2.rpt")){
    download.file("http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt", destfile="MRK_List2.rpt")
  }
  #take a subset for testing
  test_exp_set <- cortex_mrna$exp[1:10,1:50] #causes error when no mismatch
  #Add in fake gene data for which the gene could have issues with excel:
  #MARCH1 gene would go to Mar-01
  rownames(test_exp_set)[7]<-"Mar-01"
  #function should warn the user about this -if warning returned function worked
  warning_return <- 
    tryCatch(EWCE::fix.bad.mgi.symbols(test_exp_set,
                                       mrk_file_path="MRK_List2.rpt"),
             error=function(e) e, 
             warning=function(w) w)
  
  #Now test if a synonym of a gene in the list is added 
  #function should combine them and give sum reads for each sample
  # alt symbol for Tspan12 is Tm4sf12
  rownames(test_exp_set)[7]<-"Tm4sf12"
  EWCE_return = EWCE::fix.bad.mgi.symbols(test_exp_set,
                                          mrk_file_path="MRK_List2.rpt")
  #Now the returned value should be the average of the two
  sum_exp <- 
    colSums(test_exp_set[rownames(test_exp_set) %in% c("Tspan12","Tm4sf12"),])
  
  #lastly check nothing changes when there are no issues
  test_exp_set <- test_exp_set[1:5,]
  EWCE_output_same_input <- EWCE::fix.bad.mgi.symbols(test_exp_set,
                            mrk_file_path="MRK_List2.rpt")
  
  # fail if any of the 3 tests fail
  expect_equal(
    #Test 1
    all(is(warning_return,"warning"),
    #Test 2
    all(sum_exp==EWCE_return[rownames(EWCE_return)=="Tspan12",]),
    #Test 3
    all.equal(EWCE_output_same_input,test_exp_set))
    , TRUE)
})