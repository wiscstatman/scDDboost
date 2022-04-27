context("calD")
test_that("pdd gives correct result",{
    data(sim_dat)
    data_counts <- assays(sim_dat)$counts
    conditions <- colData(sim_dat)$conditions
    rownames(data_counts) <- seq_len(1000)
    bp <- BiocParallel::MulticoreParam(2)
    D_c <- calD(data_counts,bp)
    expect_equal(as.numeric(D_c[100,101]),0.01811737,tolerance=1e-6)
})
