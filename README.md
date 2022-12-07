# MCMCMixBimodal

## Introduction

MCMCMixBimodal is a package that allows the user to to generate random variables from from bimodal Densities with highly spaced modes using Metropolis-Hastings Algorithm with a specially customized proposal Density. More specifically, it allows the user to generate random variables from a mixture of two Gaussian Densities by specifying the means and variances of the component normal densities and also the mixing parameter. It also allows the user to study the convergence rate of the algorithm using relevant plots.

<!-- badges: start -->
[![R-CMD-check](https://github.com/SrijatoBhattacharyya/MCMCMixBimodal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SrijatoBhattacharyya/MCMCMixBimodal/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/SrijatoBhattacharyya/MCMCMixBimodal/branch/main/graph/badge.svg)](https://app.codecov.io/gh/SrijatoBhattacharyya/MCMCMixBimodal?branch=main)
<!-- badges: end -->


## Installation

``` r
devtools::install_github("SrijatoBhattacharyya/MCMCMixBimodal")
```

## Future work

Following this, I plan to complete the MCMCrateGaussian function. I also intend to optimize the efficiency of the functions and move the codes to c++ and finish documentation. Time permitting, I intend to modify the functions to simulate random variables from Bimodal mixture of two Cauchy Distributions by modifying the proposal density in the Metropolis-Hastings Algorithm accordingly.

## References

W.R. Gilks, S. Richardson, and D.J. Spiegelhalter, ed. (1996), Markov chain Monte Carlo in practice. Chapman and Hall, London.
