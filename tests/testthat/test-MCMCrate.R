test_that("MCMCrate works", {

  #expect warning message
  testthat::expect_error(MCMCrate(m1 = 20, m2 = 0, s1 = 1, s2 = 2, alpha = 0.5, Nsim = 10))

  #expect warning message
  testthat::expect_error(MCMCrate(m1 = 20, m2 = 0, s1 = 1, s2 = -2, alpha = 0.5, Nsim = 1000))

  #expect warning message
  testthat::expect_error(MCMCrate(m1 = 20, m2 = 0, s1 = 1, s2 = 2, alpha = 5, Nsim = 1000))

  #expect warning message
  testthat::expect_error(MCMCrate(m1 = 20, m2 = 0, s1 = 1, s2 = 2, alpha = -3, Nsim = 1000))

  #expect a vector of length 100
  testthat::expect_length(MCMCrate(m1 = 20, m2 = 0, s1 = 1, s2 = 2, alpha = 0.5, Nsim = 100), 100)


})
