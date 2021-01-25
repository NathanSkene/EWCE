# test generate.bootstrap.plots and generate.bootstrap.plots.for.transcriptome
test_that("Filter genes without 1 to 1 homolog test", {
  #use vignette data
  data("ctd")
  data("mouse_to_human_homologs")
  m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
  mouse.hits = unique(m2h[m2h$HGNC.symbol %in% example_genelist,"MGI.symbol"])
  #mouse.bg  = unique(setdiff(m2h$MGI.symbol,mouse.hits))
  mouse.bg  = unique(m2h$MGI.symbol)
  reps=100 # <- Use 100 bootstrap lists so it runs quickly, for publishable analysis use >10000
  level=1 # <- Use level 1 annotations (i.e. Interneurons)
  full_results = bootstrap.enrichment.test(sct_data=ctd,hits=mouse.hits,bg=mouse.bg,
                                           reps=reps,annotLevel=level)
  options(warn=-1)#turn off warnings for plot warning
  generate.bootstrap.plots(sct_data=ctd,hits=mouse.hits,
                           bg=mouse.bg,reps=reps,annotLevel=1,
                           full_results=full_results,listFileName="VignetteGraphs")
  options(warn=0)
  #check the BootstrapPlots folder exists and is non-empty
  test1 <- dir.exists("~/BootstrapPlots") && length(list.files("~/BootstrapPlots"))>0
  #remove folder once tested
  unlink("~/BootstrapPlots",recursive = T)
  
  data(tt_alzh)
  tt_results = ewce_expression_data(sct_data=ctd,tt=tt_alzh,annotLevel=1,
                                    ttSpecies="human",sctSpecies="mouse")
  options(warn=-1)#turn off warnings for plot warning
  full_results = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=tt_alzh,
                                                            annotLevel=1,full_results=tt_results,
                                                            listFileName="examples",
                                                            reps=reps,ttSpecies="human",
                                                            sctSpecies="mouse",onlySignif=FALSE)
  options(warn=0)
  #check the BootstrapPlots folder exists and is non-empty
  test2 <- dir.exists("BootstrapPlots") && length(list.files("BootstrapPlots"))>0
  #remove folder once tested
  unlink("BootstrapPlots",recursive = T)
  
  #fail if either function didn't create directory and add files
  expect_equal(all(test1,test2),TRUE)
})