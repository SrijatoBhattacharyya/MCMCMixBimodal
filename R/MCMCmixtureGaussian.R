#---------------------------------------------#
# function for target bimodal density
#---------------------------------------------#


target_density <- function(x, mean1, mean2, var1, var2, alpha) {
  W <- alpha * stats::dnorm(x, mean1, var1)
  V <- (1 - alpha) * stats::dnorm(x, mean2, var2)
  return(W + V)
}

#-----------------------------------------------#
# function to generate Markov chain
#-----------------------------------------------#

#' chain
#'
#' chain simulates random variables from a Markov Chain having stationery distribution as
#' a univariate Bimodal Gaussian Mixture Distribution. The means and variances of the 2
#' Gaussian distributions and the mixing parameter that form the target mixture distribution will be taken
#' as input from the user, along with the number of values to be simulated. The function returns
#' the values generated from the chain
#'
#'
#'
#'
#' @param mean1 mean of 1st univariate gaussian distribution
#' @param mean2 mean of 2nd univariate gaussian distribution
#' @param var1  variance of 1st univariate gaussian distribution
#' @param var2  variance of 2nd univariate gaussian distribution
#' @param alpha mixture proportion of the normal density with smaller mean
#' @param t starting value of the chain
#' @param Nsim  number of values to be simulated from the mixture distribution
#'
#'
#' @return vector of simulated values from the target distribution
#' @export
#'
#' @examples
#' chain(2, 1, 10, 2, 3, 0.3, 1000)
#'

chain <- function(t, mean1, mean2, var1, var2, alpha, Nsim) {
  X <- rep(0, Nsim)
  mean11 <- min(mean1, mean2)
  mean22 <- max(mean1, mean2)

  X[1] <- t # Initialized the chain

  for (i in 2:Nsim)
  {
    # Code for Proposal Density
    W <- stats::rnorm(1, (alpha * mean22) + ((1 - alpha) * X[i - 1]), 5)
    V <- stats::rnorm(1, (alpha * mean11) + ((1 - alpha) * X[i - 1]), 5)

    # Adjusting the proposal density based on the location of the previously iterated value
    if (abs(X[i - 1] - mean11) < abs(X[i - 1] - mean22)) {
      Y <- W
    } else {
      Y <- V
    }


    rho <- target_density(Y, mean1, mean2, var1, var2, alpha) / target_density(X[i - 1], mean1, mean2, var1, var2, alpha) # defined the rate parameter for Metropolis Hastings Algorithm

    # stored the next iterated value of the chain
    if (stats::runif(1) < rho) {
      X[i] <- Y
    } else {
      X[i] <- X[i - 1]
    }
  }
  return(X)
}





#--------------------------------------------#
# function to return values from simultanious running of the chain and generate histogram from the data
#--------------------------------------------#


#' MCMC mixture Gaussian
#'
#' MCMCmixtureGaussian function simulate
#' random variables from a Univariate Bimodal Gaussian Mixture Distribution. The means and variances of the 2
#' Gaussian distributions and the mixing parameter that form the target mixture distribution will be taken
#' as input from the user, along with the number of values to be simulated. This function generates 1000 replications
#' of the markov chain with stationery distribution as the target distribution and plots the histogram of the values produced.
#'
#' @param mean1 mean of 1st univariate gaussian distribution
#' @param mean2 mean of 2nd univariate gaussian distribution
#' @param var1  variance of 1st univariate gaussian distribution
#' @param var2  variance of 2nd univariate gaussian distribution
#' @param alpha mixture proportion of the normal density with smaller mean
#' @param Nsim  number of values to be simulated from the mixture distribution
#'
#' @return vector of simulated values from the target distribution
#' @export
#'
#' @examples
#' MCMCmixtureGaussian(1, 10, 2, 3, 0.3, 1000)
#'
MCMCmixtureGaussian <- function(mean1, mean2, var1, var2, alpha, Nsim) {
  mean11 <- min(mean1, mean2)
  mean22 <- max(mean1, mean2)

  n <- 1000 # declared number of replications of the chain


  Z <- do.call("rbind", replicate(n, chain(mean11, mean1, mean2, var1, var2, alpha, Nsim), simplify = F))
  R <- Z[, 500]

  # histogram of generated values along with the density curve
  graphics::hist(R, breaks = 70, prob = T)
}
