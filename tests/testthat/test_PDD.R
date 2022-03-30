context("PDD")
test_that("PDD gives correct result",{
    data(sim_dat)
    data_counts <- assays(sim_dat)$counts
    conditions <- colData(sim_dat)$conditions
    rownames(data_counts) <- seq_len(1000)
    D_c <- cal_D(data_counts,1)
    pdd <- PDD(data = data_counts,cd = conditions, ncores = 1, D = D_c)
    expect_equal(pdd[101],0.005305089,tolerance=1e-6)
})
