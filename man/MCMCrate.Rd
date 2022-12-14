% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMCrate.R
\name{MCMCrate}
\alias{MCMCrate}
\title{MCMC rate}
\usage{
MCMCrate(m1, m2, s1, s2, alpha, Nsim, density = "Normal")
}
\arguments{
\item{m1}{(numeric data type): mean of 1st univariate gaussian distribution or median of 1st cauchy Distribution depending on value of density parameter}

\item{m2}{(numeric data type) : mean of 2nd univariate gaussian distribution or median of 2nd cauchy Distribution depending on value of density parameter}

\item{s1}{(numeric data type) : variance of 1st univariate gaussian distribution or scale parameter of 1st Cauchy Distribution depending on value of density parameter}

\item{s2}{(numeric data type) : variance of 2nd univariate gaussian distribution or scale parameter of 2nd Cauchy Distribution depending on value of density parameter}

\item{alpha}{(numeric data type) : value of mixing parameter}

\item{Nsim}{number of simulations of Markov Chain having stationery distribution as the target mixture distribution}

\item{density}{(character data type) : specifies the density function as either Normal or Cauchy. Default = "Normal"}
}
\value{
Nsim-long vector of supremum of the absolute difference between the empirical
CDF of the 1000 generated values of the chain and the actual CDF of the
target distribution at each iteration.
}
\description{
MCMCrate generates 1000 Markov Chain values simultaneously
at each iteration of the Algorithm. These chains have stationery distribution
same as the target mixture distribution. This function generates Nsim values for each chain. As its output, it will return the Nsim-long vector of supremum of the absolute
difference between the empirical CDF of the 1000 generated values of the chain and the actual CDF of the
target distribution at each iteration. It will also plot these differences for the first 50 iterations.
}
\examples{
MCMCrate(m1 = 20, m2 = 0, s1 = 1, s2 = 2, alpha = 0.5, Nsim = 1000)
MCMCrate(m1 = 10, m2 = 0, s1 = 2, s2 = 4, alpha = 0.5, Nsim = 1000, density = "Cauchy")

}
\references{
W.R. Gilks, S. Richardson, and D.J. Spiegelhalter, ed. (1996),
Markov chain Monte Carlo in practice. Chapman and Hall, London.
}
