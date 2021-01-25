#Test for add.res.to.merging.list, using sample data in vignette
test_that("merging EWCE results from multiple sets", {
  #load vignette data
  # Load some vignette data
  data(tt_alzh)
  data(tt_alzh_BA36)
  data(tt_alzh_BA44)
  data("ctd")
  
  # Run EWCE analysis
  tt_results = ewce_expression_data(sct_data=ctd,tt=tt_alzh,annotLevel=1,
                                    ttSpecies="human",sctSpecies="mouse")
  tt_results_36 = ewce_expression_data(sct_data=ctd,tt=tt_alzh_BA36,
                                       annotLevel=1,ttSpecies="human",
                                       sctSpecies="mouse")
  tt_results_44 = ewce_expression_data(sct_data=ctd,tt=tt_alzh_BA44,
                                       annotLevel=1,ttSpecies="human",
                                       sctSpecies="mouse")
  # Fill a list with the results
  results = add.res.to.merging.list(tt_results)
  results = add.res.to.merging.list(tt_results_36,results)
  results = add.res.to.merging.list(tt_results_44,results)
  
  #check returns
  test1 <- length(results)==6
  #make sure all metrics returned for each
  test2 <- all.equal(unlist(lapply(results, function(x) names(x))),
                 rep(c("Direction","bootstrap_data","hitCells"),6))
  #check output equal to input before combination
  test3 <- all.equal(results[[1]]$hitCells,tt_results$hit.cells.up)
  test4 <- all.equal(results[[5]]$hitCells,tt_results_44$hit.cells.up)
  test5 <- all.equal(tt_results_36$bootstrap_data.up,results[[3]]$bootstrap_data)
  
  # Perform the merged analysis
  merged_res = merged_ewce(results,reps=10) # <- For publication reps should be higher
  
  ewce_plot_res <- ewce.plot(merged_res)$plain
  # fail if any but ggplot returned
  test6 <- all(class(ewce_plot_res)==c("gg","ggplot"))
  
  # fail if any subtest isn't true
  expect_equal(all(test1,test2,test3,test4,test5,test6),TRUE)
})
