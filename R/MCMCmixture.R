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
#' target_density(x = 7, m1 = 5, m2 = 1, s1 = 3, s2 = 2, alpha = 0.3, density = "Cauchy")
#'
target_density <- function(x, m1, m2, s1, s2, alpha, density = "Normal") {

  #check that scale parameter of 1st distribution is positive
  if(s1 <= 0 ){
    stop("Scale parameter cannot be non-positive; supply new scale parameter")
  }

  #check that scale parameter of 2nd distribution is positive
  if(s2 <0 ){
    stop("Scale parameter cannot be non-positive; supply new scale parameter")
  }

  #check that the mixing parameter is not greater than 1
  if(alpha > 1 ){
    stop("mixing parameter cannot be greater than 1")
  }

  #check that the mixing parameter is not negative
  if(alpha <0 ){
    stop("mixing parameter cannot be negative")
  }

  #Calculation for density for mixture of Cauchy distributions
  if (density == "Cauchy") {
    W <- alpha * stats::dcauchy(x, m1, s1)
    V <- (1 - alpha) * stats::dcauchy(x, m2, s2)
  } else {

    #Calculation for density for mixture of Normal distributions
    W <- alpha * stats::dnorm(x, m1, s1)
    V <- (1 - alpha) * stats::dnorm(x, m2, s2)
  }

  #retrun the value of density
  return(W + V)
}






#-----------------------------------------------#
# function to generate Markov chain
#-----------------------------------------------#

#' chain
#'
#' chain function simulates random variables from a Markov Chain having stationery distribution as
#' a univariate Bimodal Gaussian/Cauchy Mixture Distribution. The means/medians and variances/scale parameters of the 2
#' Gaussian/Cauchy distributions and the mixing parameter that form the target mixture distribution will be taken
#' as input from the user, along with the number of values to be simulated. The function returns
#' the values generated from the chain
#'
#'
#'
#' @param t starting value of the chain
#' @inheritParams target_density
#' @param Nsim  number of simulations of Markov Chain having stationery distribution as the target mixture distribution
#'
#'
#' @return vector of simulated values from the target distribution
#' @export
#'
#' @references
#'
#' W.R. Gilks, S. Richardson, and D.J. Spiegelhalter, ed. (1996),
#' Markov chain Monte Carlo in practice. Chapman and Hall, London.
#'
#' @examples
#' chain(2, 1, 10, 2, 3, 0.3, 1000)
#' chain(5, 4, 3, 3, 3, 0.5, 1000)
#' chain(-3, -2, -4, 2, 2, 0.9, 1000, density = "Cauchy")
#'
chain <- function(t, m1, m2, s1, s2, alpha, Nsim, density = "Normal") {

  #Check that Nsim is greater than 50 and that Nsim is a natural number
  if(Nsim <= 50 || Nsim %% 1 != 0){
    stop("Please resupply positive integer value of Nsim that is greater than 50 ")
  }

  #Defined vector to store the chain values
  X <- rep(0, Nsim)

  #calculate the minimum of the 2 input location parameters
  m11 <- min(m1, m2)
  m22 <- max(m1, m2)

  #identifying scale parameter of distributions with lower and higher location paramters
  if (m11 == m1) {
    s11 <- s1
    s22 <- s2
  } else {
    s11 <- s2
    s22 <- s1
  }

  X[1] <- t # Initialized the chain

  #Running the Chain
  for (i in 2:Nsim)
  {
    # Code for Proposal Density
    if (density == "Cauchy") {
      W <- stats::rcauchy(1, (alpha * m22) + ((1 - alpha) * X[i - 1]), s2 + 10)
      V <- stats::rcauchy(1, (alpha * m11) + ((1 - alpha) * X[i - 1]), s1 + 10)
    } else {
      W <- stats::rnorm(1, (alpha * m22) + ((1 - alpha) * X[i - 1]), s2 + 10)
      V <- stats::rnorm(1, (alpha * m11) + ((1 - alpha) * X[i - 1]), s1 + 10)
    }

    # Adjusting the proposal density based on the location of the previously iterated value
    if (abs(X[i - 1] - m11) < abs(X[i - 1] - m22)) {
      Y <- W
    } else {
      Y <- V
    }

    #calculating acceptance probability for Metropolis-Hastings Algo
    rho <- target_density(Y, m1, m2, s1, s2, alpha) / target_density(X[i - 1], m1, m2, s1, s2, alpha) # defined the rate parameter for Metropolis Hastings Algorithm

    # stored the next iterated value of the chain
    if (stats::runif(1) < rho) {
      X[i] <- Y
    } else {
      X[i] <- X[i - 1]
    }
  }

  #returned the vector of simulated values
  return(X)
}





#--------------------------------------------#
# function to return values from simultanious running of the chain and generate histogram from the data
#--------------------------------------------#


#' MCMC mixture
#'
#' MCMCmixture function simulates
#' random variables from a Univariate Bimodal Gaussian/Cauchy Mixture Distribution. The means/medians and variances/scale parameters of the 2
#' Gaussian/Cauchy distributions and the mixing parameter that form the target mixture distribution will be taken
#' as input from the user, along with the iteration number of the Markov Chains to be simulated. This function generates 1000 replications
#' of the markov chain with stationery distribution as the target distribution, returns the 1000-long vector of generated values at the Nsim-th
#' iteration of each chain and plots the histogram of the values produced.
#'
#'
#' @inheritParams target_density
#' @param Nsim  number of simulations of Markov Chain having stationery distribution as the target mixture distribution
#'
#' @return vector of simulated values from the target distribution
#' @export
#'
#'#' @references
#'
#' W.R. Gilks, S. Richardson, and D.J. Spiegelhalter, ed. (1996),
#' Markov chain Monte Carlo in practice. Chapman and Hall, London.
#'
#' @examples
#' MCMCmixture(1, 10, 2, 3, 0.3, 1000, density = "Normal")
#' MCMCmixture(-10, 20, 10, 5, 0.6, 2000, density = "Normal")
#' MCMCmixture(-1, 7, 1, 3, 0.5, 3000, density = "Cauchy")
#'
MCMCmixture <- function(m1, m2, s1, s2, alpha, Nsim, density = "Normal") {

  #calculate the minimum of the 2 input location parameters
  m11 <- min(m1, m2)
  m22 <- max(m1, m2)


  # declared number of replications of the chain
  n <- 1000

  #replicating the chain n times
  Z <- do.call("rbind", replicate(n, chain(m11, m1, m2, s1, s2, alpha, Nsim, density), simplify = F))

  # Collecting the Nsimth values in each chain
  Value <- Z[, Nsim]
  x <- min(Value) - 3:max(Value) + 3

  # histogram of generated values along with the density curve
  graphics::hist(Value, breaks = 70, prob = T, main = paste("Histogram of simulated values") )
  graphics::curve(target_density(x, m1, m2, s1, s2, alpha, density), min(Value), max(Value), add = TRUE)

  #returning collection of the Nsimth values in each chain
  return(Value)
}
