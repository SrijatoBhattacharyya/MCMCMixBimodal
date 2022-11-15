#' MCMC rate Gaussian
#'
#'
#' MCMCrateGaussian generate 500 Markov Chain values simultaneously
#' at each iteration of the Algorithm. As its output, it will return the curve of supremum of the absolute
#' difference between the empirical CDF of the 500 generated values of the chain and the actual CDF of the
#' target distribution at each iteration.
#'
#' @param mean1 mean of 1st univariate gaussian distribution
#' @param mean2 mean of 2nd univariate gaussian distribution
#' @param var1  variance of 1st univariate gaussian distribution
#' @param var2  variance of 2nd univariate gaussian distribution
#' @param nsim  number of values to be simulated from the mixture distribution
#'
#' @return plot of supremum of the absolute
#' difference between the empirical CDF of the 500 generated values of the chain and the actual CDF of the
#' target distribution at each iteration.
#' @export
#'
#' @examples MCMCrateGaussian(mean1, mean2, var1, var2, nsim)
#'
#'
MCMCrateGaussian <- function(mean1, mean2, var1, var2, nsim){

  return(0)

}
