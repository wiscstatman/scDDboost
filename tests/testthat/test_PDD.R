context("pdd")
test_that("pdd gives correct result",{
    data(sim_dat)
    data_counts <- assays(sim_dat)$counts
    conditions <- colData(sim_dat)$conditions
    rownames(data_counts) <- seq_len(1000)
    bp <- BiocParallel::MulticoreParam(2)
    D_c <- calD(data_counts,bp)
    pdd <- pdd(data = data_counts,cd = conditions, bp, D = D_c, random = FALSE)
    expect_equal(as.numeric(pdd[101]),0.008143924,tolerance=1e-6)
})
