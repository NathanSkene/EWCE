# Test for mean calculation from calculate_meanexp_for_level in generate_celltype_data
test_that("SingleCellExperiment data type matches original input format results", {
  #SCT dataset
  cortex_mrna <- ewceData::cortex_mrna()
  cortex_mrna$exp <- cortex_mrna$exp[1:6000,]#reduce number of genes for speed
  #Make SCE object from SCT
  cortex_mrna_SCE <-
    SingleCellExperiment::SingleCellExperiment(
      assays=list(counts=cortex_mrna$exp),
      colData=cortex_mrna$annot)

  #Ensure output is the same for the two data types
  cortex_mrna$exp = suppressWarnings(fix_bad_mgi_symbols(cortex_mrna$exp))
  cortex_mrna_SCE = suppressWarnings(fix_bad_mgi_symbols(cortex_mrna_SCE))
  #fail if not the exact same
  expect_equal(all.equal(assays(cortex_mrna_SCE)$counts,cortex_mrna$exp),TRUE)
})
