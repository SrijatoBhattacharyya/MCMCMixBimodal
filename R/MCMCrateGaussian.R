
# source(here::here("R", "MCMCmixtureGaussian.R"))



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
#' @param alpha value of mixing parameter
#' @param Nsim  number of values to be simulated from the mixture distribution
#'
#' @return plot of supremum of the absolute
#' difference between the empirical CDF of the 500 generated values of the chain and the actual CDF of the
#' target distribution at each iteration.
#' @export
#'
#' @examples
#' MCMCrateGaussian(mean1 = 20, mean2 = 0, var1= 1, var2 = 2, alpha= 0.5, Nsim= 1000)
#'



MCMCrateGaussian <- function(mean1, mean2, var1, var2, alpha, Nsim) {
  n = 1000
  S=c()

  Z <- do.call("rbind", replicate(n, chain(mean1, mean1, mean2, var1, var2, alpha, Nsim), simplify = F))

  # Z=do.call("rbind",replicate(n, chain( t = mean1, mean1 = mean1, mean2 = mean2, var1 = var1, var2 = var2, alpha = alpha, Nsim = Nsim), simplify = F))
  for(j in 1: 50)
  {
    P=Z[,j]
    V=sort(P)
    E=c()
    for(i in 1:n)
    {
      E[i]=i/n
    }
    O=c()
    for(i in 1:n)
    {
      x = V[i]
      O[i]=abs(E[i]- ((alpha * stats::pnorm(x, mean1, var1))+((1 - alpha) * stats::pnorm(x, mean2, var2))) )
    }
    M=max(O)
    S[j]=M
  }

  plot(S)

}

