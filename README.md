# MCMCMixBimodal

## Introduction

MCMCMixBimodal is a package that allows the user to to generate random variables from from bimodal Densities with highly spaced modes using Metropolis-hastings Algorithm with a specially customized proposal Density. More specifically, it allows the user to generate random variables from a mixture of two Gaussian Densities by specifying the means and variances of the component normal densities and also the mixing parameter. It also allows the user to study the convergence rate of the algorithm using relevant plots.

## Installation

``` r
devtools::install_github("SrijatoBhattacharyya/MCMCMixBimodal")
```

## Future work

Following this, I plan to complete the MCMCrateGaussian function. I also intend to optimize the efficiency of the functions and move the codes to c++ and finish documentation. Time permitting, I intend to modify the functions to simulate random variables from Bimodal mixture of two Cauchy Distributions by modifying the proposal density in the Metropolis-Hastings Algorithm accordingly.
