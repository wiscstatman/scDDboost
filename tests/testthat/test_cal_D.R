context("cal_D")
test_that("PDD gives correct result",{
    data(sim_dat)
    data_counts <- assays(sim_dat)$counts
    conditions <- colData(sim_dat)$conditions
    rownames(data_counts) <- seq_len(1000)
    D_c <- cal_D(data_counts,1)
    expect_equal(as.numeric(D_c[100,101]),0.01811737,tolerance=1e-6)
})
