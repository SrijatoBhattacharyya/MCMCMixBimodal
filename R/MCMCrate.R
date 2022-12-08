
# source(here::here("R", "MCMCmixture.R"))



#' MCMC rate
#'
#'
#' MCMCrate generate 1000 Markov Chain values simultaneously
#' at each iteration of the Algorithm and generates Nsim values for each chain. As its output, it will return the Nsim-long vector of supremum of the absolute
#' difference between the empirical CDF of the 1000 generated values of the chain and the actual CDF of the
#' target distribution at each iteration. It will also plot these differences.
#'
#' @inheritParams target_density
#' @param Nsim  number of values to be simulated from the mixture distribution
#'
#' @return vector of supremum of the absolute difference between the empirical
#' CDF of the 500 generated values of the chain and the actual CDF of the
#' target distribution at each iteration.
#' @export
#'
#' @examples
#' MCMCrate(m1 = 20, m2 = 0, s1 = 1, s2 = 2, alpha = 0.5, Nsim = 1000)
#' MCMCrate(m1 = 10, m2 = 0, s1 = 2, s2 = 4, alpha = 0.5, Nsim = 1000, density = "Cauchy")
#'
MCMCrate <- function(m1, m2, s1, s2, alpha, Nsim, density = "Normal") {
  n <- 1000
  S <- c()

  Z <- do.call("rbind", replicate(n, chain(m1, m1, m2, s1, s2, alpha, Nsim, density), simplify = F))

  for (j in 1:Nsim)
  {
    P <- Z[, j]
    V <- sort(P)
    E <- c()
    for (i in 1:n){
      E[i] <- i / n
    }
    O <- c()
    if (density == "Cauchy") {
      for (i in 1:n){
        x <- V[i]
        O[i] <- abs(E[i] - ((alpha * stats::pcauchy(x, m1, s1)) + ((1 - alpha) * stats::pcauchy(x, m2, s2))))
      }
    }
    else {
      for (i in 1:n)
      {
        x <- V[i]
        O[i] <- abs(E[i] - ((alpha * stats::pnorm(x, m1, s1)) + ((1 - alpha) * stats::pnorm(x, m2, s2))))
      }
    }

    M <- max(O)
    S[j] <- M
  }

  Supnorm = S[1:50]
  plot(Supnorm)
  return(S)
}
