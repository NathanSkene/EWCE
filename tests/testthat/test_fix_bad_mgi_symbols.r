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
  
  #catch error when no exp inputted
  error_return <- 
    tryCatch(EWCE::fix.bad.mgi.symbols(exp=NULL,
                                       mrk_file_path="MRK_List2.rpt"),
             error=function(e) e, 
             warning=function(w) w)
  
  #catch error when no test_exp_set inputted as character
  test_exp_set_char <- apply(test_exp_set, 2, as.character)
  rownames(test_exp_set_char) <- rownames(test_exp_set)
  error_return2 <- 
    tryCatch(EWCE::fix.bad.mgi.symbols(exp=test_exp_set_char,
                                       mrk_file_path="MRK_List2.rpt"),
             error=function(e) e, 
             warning=function(w) w)
  
  #pass data table rather than data frame or matrix
  error_return3 <- 
    tryCatch(EWCE::fix.bad.mgi.symbols(exp=data.table::data.table(as.data.frame(test_exp_set)),
                                       mrk_file_path="MRK_List2.rpt"),
             error=function(e) e, 
             warning=function(w) w)
  
  #function should warn the user about this -if warning returned function worked
  warning_return <- 
    tryCatch(EWCE::fix.bad.mgi.symbols(test_exp_set,
                                       mrk_file_path="MRK_List2.rpt"),
             error=function(e) e, 
             warning=function(w) w)
  
  #running on hgnc rather than mgi should return warnings
  warning_return2 <- 
    tryCatch(EWCE::fix.bad.hgnc.symbols(test_exp_set),
             error=function(e) e, 
             warning=function(w) w)
  
  #running on hgnc rather than mgi should return warnings
  warning_return3 <- 
    tryCatch(EWCE::fix.bad.hgnc.symbols(data.table::data.table(as.data.frame(test_exp_set))),
             error=function(e) e, 
             warning=function(w) w)
  
  options(warn=-1)
  hgnc_return <- EWCE::fix.bad.hgnc.symbols(test_exp_set)
  options(warn=0)
  
  #Now test if a synonym of a gene in the list is added 
  #function should combine them and give sum reads for each sample
  # alt symbol for Tspan12 is Tm4sf12
  rownames(test_exp_set)[7]<-"Tm4sf12"
  EWCE_return = EWCE::fix.bad.mgi.symbols(test_exp_set,
                                          mrk_file_path="MRK_List2.rpt",
                                          printAllBadSymbols=T)
  #Now the returned value should be the average of the two
  sum_exp <- 
    colSums(test_exp_set[rownames(test_exp_set) %in% c("Tspan12","Tm4sf12"),])
  
  #check nothing changes when there are no issues
  test_exp_set <- test_exp_set[1:5,]
  EWCE_output_same_input <- EWCE::fix.bad.mgi.symbols(test_exp_set,
                            mrk_file_path="MRK_List2.rpt")
  
  
  #check input runs on large size of incorrect mgi gene names
  #check it catches warning still but doesn't give an error
  warning_return4 <- 
    tryCatch(EWCE_output_large <- EWCE::fix.bad.mgi.symbols(cortex_mrna$exp[1:10000,1:1000],
                                                            mrk_file_path="MRK_List2.rpt"),
             error=function(e) e, 
             warning=function(w) w)
  
  #remove folder once tested
  unlink("MRK_List2.rpt",recursive = T)
  
  # fail if any of the 3 tests fail
  expect_equal(
    #Test 0.1
    all(is(error_return,"error"),
    #Test 0.2
    is(error_return2,"warning"),
    #Test 0.3
    is(error_return3,"error"),
    #Test 1
    is(warning_return,"warning"),
    #Test 1.2
    is(warning_return2,"warning"),
    #Test 1.3
    is(warning_return3,"error"),
    #Test 1.4
    class(hgnc_return)[1]=="matrix",
    #Test 2
    all(sum_exp==EWCE_return[rownames(EWCE_return)=="Tspan12",]),
    #Test 3
    all.equal(EWCE_output_same_input,test_exp_set),
    #Test 4
    is(warning_return4,"warning"))
    , TRUE)
})