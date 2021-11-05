test_that("DGE works", {
    
    if(!is_32bit()){
        set.seed(1234)
        cortex_mrna <- ewceData::cortex_mrna()
        n_genes <- 100
        exp <- cortex_mrna$exp[seq(1,n_genes),]
        level2annot <- cortex_mrna$annot$level2class
        
        #### limma ####
        limma_res <- EWCE::: run_limma(exp = exp,
                                      level2annot = level2annot)
        testthat::expect_true(
            all(c("coefficients","rank","assign","F","p.value") %in% names(limma_res)))
        testthat::expect_true(EWCE:::is_matrix(limma_res$coefficients))
        testthat::expect_equal(sum(limma_res$q<.05),100)
        
        #### DESeq2 ####
        deqseq2_res <- EWCE:::run_deseq2(exp = exp,
                                         level2annot =  level2annot,
                                         test = "LRT")
        testthat::expect_true(
            all(c("baseMean","log2FoldChange","pvalue","padj") %in% names(deqseq2_res)))
        testthat::expect_length(deqseq2_res$padj, n_genes)
        testthat::expect_equal(nrow(subset(deqseq2_res, padj<.05)),80)
        
        #### MAST ####
        mast_res <- EWCE:::run_mast(exp = exp,
                                    level2annot = level2annot)
        testthat::expect_true(
            all(c("Class","lrstat","p.value") %in% names(mast_res)))
        testthat::expect_true(is.data.frame(mast_res))
        testthat::expect_equal(nrow(subset(mast_res, q<.05)),1917)
    }
})
