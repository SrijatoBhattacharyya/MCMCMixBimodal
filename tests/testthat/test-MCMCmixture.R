testthat::test_that("MCMCmixture test cases", {

  #expecting a value
  testthat::expect_length(target_density(x = 8, m1 = 2, m2 = 7, s1 = 2, s2 = 1, alpha = 0.1, density = "Normal"),1)

  #expecting warning message
  testthat::expect_error(target_density(x = 8, m1 = 2, m2 = 7, s1 = -2, s2 = 1, alpha = 0.1, density = "Normal"))

  #expecting warning message
  testthat::expect_error(target_density(x = 8, m1 = 2, m2 = 7, s1 = 2, s2 = 1, alpha = 10, density = "Normal"))

  #expecting warning message
  testthat::expect_error(target_density(x = 8, m1 = 2, m2 = 7, s1 = 2, s2 = 1, alpha = -5, density = "Normal"))


  # expecting a vector of length 100
  testthat::expect_length(chain(2, 1, 10, 2, 3, 0.3, 100), 100)

  # expecting a vector of length 1000
  testthat::expect_length(chain(2, 1, 10, 2, 3, 0.3, 1000), 1000)

  #expecting warning message
  testthat::expect_error(chain(2, 1, 10, 2, 3, 3, 1000))

  #expecting warning message
  testthat::expect_error(chain(2, 1, 10, -2, 3, 0.3, 1000))


  # expecting warning message
  testthat::expect_error(MCMCmixture(1, 10, 2, 3, 3, -1000, density = "Normal"))



  # expecting a vector of length 1000
  testthat::expect_length(MCMCmixture(1, 10, 2, 3, 0.3, 1000, density = "Normal"), 1000)

  # expecting a vector of length 1000
  testthat::expect_length(MCMCmixture(1, 10, 2, 2, 0.5, 1000, density = "Cauchy"), 1000)


  # expecting warning message
  testthat::expect_error(MCMCmixture(1, 10, 2, 3, 3, 1000, density = "Normal"))

  # expecting warning message
  testthat::expect_error(MCMCmixture(1, 10, 2, 3, 0.3, 10, density = "Normal"))

  # expecting warning message
  testthat::expect_error(MCMCmixture(1, 10, -2, 3, 3, 1000, density = "Normal"))

  # expecting warning message
  testthat::expect_error(MCMCmixture(1, 10, 2, 3, 3, 100.5, density = "Normal"))




})
