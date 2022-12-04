#---------------------------------------------#
# function for target bimodal density
#---------------------------------------------#


#' target_density
#'
#'
#' target_density gives the value of the density of a univariate Bimodal Gaussian / Bimodal Cauchy Mixture Distribution. The means/medians and variances/scale parameters of the 2
#' Gaussian/Cauchy distributions and the mixing parameter are taken as input from the user.
#'
#'
#' @param x     argument value of the density function
#' @param m1 mean of 1st univariate gaussian distribution or median of 1st cauchy Distribution depending on value of density parameter
#' @param m2 mean of 2nd univariate gaussian distribution or median of 2nd cauchy Distribution depending on value of density parameter
#' @param s1  variance of 1st univariate gaussian distribution or scale parameter of 1st Cauchy Distribution depending on value of density parameter
#' @param s2  variance of 2nd univariate gaussian distribution or scale parameter of 2nd Cauchy Distribution depending on value of density parameter
#' @param alpha value of mixing parameter
#' @param density specifies the density function as either Normal or Cauchy
#'
#' @return value of the density function at x
#' @export
#'
#'
#' @examples
#' target_density(x = 1, m1 = 5, m2 = 0, s1 = 3, s2 = 1, alpha = 0.5)
#' target_density(x = 8, m1 = 2, m2 = 7, s1 = 2, s2 = 1, alpha = 0.1, density = "Normal")
#'target_density(x = 7, m1 = 5, m2 = 1, s1 = 3, s2 = 2, alpha = 0.3, density = "Cauchy")
#'
target_density <- function(x, m1, m2, s1, s2, alpha, density = "Normal") {

  if(density == "Cauchy"){
    W <- alpha * stats::dcauchy(x, m1, s1)
    V <- (1 - alpha) * stats::dcauchy(x, m2, s2)
  }
  else{
    W <- alpha * stats::dnorm(x, m1, s1)
    V <- (1 - alpha) * stats::dnorm(x, m2, s2)
  }
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
#' @param t starting value of the chain
#' @inheritParams target_density
#' @param Nsim  number of values to be simulated from the mixture distribution
#'
#'
#' @return vector of simulated values from the target distribution
#' @export
#'
#' @examples
#' chain(2, 1, 10, 2, 3, 0.3, 1000)
#' chain(5, 4, 3, 3, 3, 0.5, 1000)
#' chain(-3, -2, -4, 2, 2, 0.9, 1000, density = "Cauchy")

chain <- function(t, m1, m2, s1, s2, alpha, Nsim, density = "Normal") {
  X <- rep(0, Nsim)
  m11 <- min(m1, m2)
  m22 <- max(m1, m2)
  if(m11 == m1){
    s11 <- s1
    s22 <- s2
  }else{
      s11 <- s2
      s22 <- s1
    }

  X[1] <- t # Initialized the chain


  for (i in 2:Nsim)
  {
    # Code for Proposal Density
    if(density == "Cauchy"){
      W <- stats::rcauchy(1, (alpha * m22) + ((1 - alpha) * X[i - 1]), s2 +  10)
      V <- stats::rcauchy(1, (alpha * m11) + ((1 - alpha) * X[i - 1]), s1 + 10)
    }
    else{
      W <- stats::rnorm(1, (alpha * m22) + ((1 - alpha) * X[i - 1]), s2 +  10)
      V <- stats::rnorm(1, (alpha * m11) + ((1 - alpha) * X[i - 1]), s1 + 10)
    }

    # Adjusting the proposal density based on the location of the previously iterated value
    if (abs(X[i - 1] - m11) < abs(X[i - 1] - m22)) {
      Y <- W
    } else {
      Y <- V
    }


    rho <- target_density(Y, m1, m2, s1, s2, alpha) / target_density(X[i - 1], m1, m2, s1, s2, alpha) # defined the rate parameter for Metropolis Hastings Algorithm

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


#' MCMC mixture
#'
#' MCMCmixture function simulate
#' random variables from a Univariate Bimodal Gaussian Mixture Distribution. The means and variances of the 2
#' Gaussian distributions and the mixing parameter that form the target mixture distribution will be taken
#' as input from the user, along with the number of values to be simulated. This function generates 1000 replications
#' of the markov chain with stationery distribution as the target distribution and plots the histogram of the values produced.
#'
#'
#' @inheritParams target_density
#' @param Nsim  number of values to be simulated from the mixture distribution
#'
#' @return vector of simulated values from the target distribution
#' @export
#'
#' @examples
#' MCMCmixture(1, 10, 2, 3, 0.3, 1000, density = "Normal")
#' MCMCmixture(-10, 20, 10, 5, 0.6, 1000, density = "Normal")
#' MCMCmixture(-10, 20, 10, 10, 0.5, 1000, density = "Cauchy")
#'
#'
MCMCmixture <- function(m1, m2, s1, s2, alpha, Nsim,  density = "Normal") {
  m11 <- min(m1, m2)
  m22 <- max(m1, m2)

  n <- 1000 # declared number of replications of the chain


  Z <- do.call("rbind", replicate(n, chain(m11, m1, m2, s1, s2, alpha, Nsim, density), simplify = F))
  R <- Z[, 500]
  x = min(R)-3:max(R)+3

  # histogram of generated values along with the density curve
  graphics::hist(R, breaks = 70, prob = T)
  graphics::curve(target_density(x, m1, m2, s1, s2, alpha,  density), min(R), max(R), add = TRUE)
  return(R)
}

