
# source(here::here("R", "MCMCmixture.R"))



#' MCMC rate
#'
#'
#' MCMCrate generates 1000 Markov Chain values simultaneously
#' at each iteration of the Algorithm. These chains have stationery distribution
#' same as the target mixture distribution. This function generates Nsim values for each chain. As its output, it will return the Nsim-long vector of supremum of the absolute
#' difference between the empirical CDF of the 1000 generated values of the chain and the actual CDF of the
#' target distribution at each iteration. It will also plot these differences for the first 50 iterations.
#'
#' @inheritParams target_density
#' @param Nsim  number of simulations of Markov Chain having stationery distribution as the target mixture distribution
#'
#' @return Nsim-long vector of supremum of the absolute difference between the empirical
#' CDF of the 1000 generated values of the chain and the actual CDF of the
#' target distribution at each iteration.
#' @export
#'
#'#' @references
#'
#' W.R. Gilks, S. Richardson, and D.J. Spiegelhalter, ed. (1996),
#' Markov chain Monte Carlo in practice. Chapman and Hall, London.
#'
#'
#' @examples
#' MCMCrate(m1 = 20, m2 = 0, s1 = 1, s2 = 2, alpha = 0.5, Nsim = 1000)
#' MCMCrate(m1 = 10, m2 = 0, s1 = 2, s2 = 4, alpha = 0.5, Nsim = 1000, density = "Cauchy")
#'
MCMCrate <- function(m1, m2, s1, s2, alpha, Nsim, density = "Normal") {

  #defined number of replications of the chain
  n <- 1000

  #defined vector to store supnorm values
  S <- c()

  #replicated the chain n times
  Z <- do.call("rbind", replicate(n, chain(m1, m1, m2, s1, s2, alpha, Nsim, density), simplify = F))


  #for loop to find supnorm values
  for (j in 1:Nsim)
  {
    P <- Z[, j]
    V <- sort(P)
    E <- c()

    #finding empirical cdf
    for (i in 1:n){
      E[i] <- i / n
    }
    O <- c()

    # calculating supnorm for Cauchy case
    if (density == "Cauchy") {
      for (i in 1:n){
        x <- V[i]
        O[i] <- abs(E[i] - ((alpha * stats::pcauchy(x, m1, s1)) + ((1 - alpha) * stats::pcauchy(x, m2, s2))))
      }
    }else {

      # calculating supnorm for Normal case
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

  #Poltting supnorm curve from 1st to 50th iterating
  plot(Supnorm)

  #returning all the supnorm values
  return(S)
}
