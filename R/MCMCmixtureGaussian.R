#' MCMC mixture Gaussian
#'
#' MCMCmixtureGaussian function simulate
#'random variables from a Univariate Gaussian Mixture Distribution. The means and variances of the 2
#'Gaussian distributions and the mixing parameter that form the target mixture distribution will be taken
#'as input from the user, along with the number of values to be simulated.
#'
#' @param mean1 mean of 1st univariate gaussian distribution
#' @param mean2 mean of 2nd univariate gaussian distribution
#' @param var1  variance of 1st univariate gaussian distribution
#' @param var2  variance of 2nd univariate gaussian distribution
#' @param nsim  number of values to be simulated from the mixture distribution
#'
#' @return vector of simulated values from the target distribution
#' @export
#'
#' @examples
#' MCMCmixtureGaussian(mean1, mean2, var1, var2, alpha, nsim)
#'
MCMCmixtureGaussian <- function(mean1, mean2, var1, var2, alpha, nsim){

  #made the function for bimodal density
  f <- function(x)
  {
    W= alpha * dnorm(x, mean1, var1)
    V=(1- alpha)*dnorm(x, mean2 ,var2)
    return(W+V)
  }
  # curve(f(x), -10, 30)
  # Nsim=500
  n=1000                              # declared number of replications of the chain
  X=rep(0, times=Nsim)                #initialized vector to store chain values

  #function to generate Markov chain

  h <- function(t)
  {
    X[1]=t
    for(i in 2:Nsim)
    {
      W=rnorm(1,( alpha * mean2) + ((1 - alpha) * X[i-1]), 5)
      V=rnorm(1,(alpha * mean1) + ((1 - alpha) * X[i-1]), 5)
      if(abs(X[i-1] - mean1) < abs(X[i-1] - mean2))
      {
        Y=W
      }
      else
      {
        Y=V
      }
      rho=f(Y) / f(X[i-1])
      if(runif(1) < rho)
      {
        X[i] = Y
      }
      else
      {
        X[i] = X[i-1]
      }
    }
    return(X)
  }
  Z=do.call("rbind",replicate(n, h(mean1), simplify = F))
  R=Z[,500]

  hist(R, breaks=70, prob=T)
  curve(f(x), mean1 - 3*var1, mean2 + 3*var2, add=T)
}

