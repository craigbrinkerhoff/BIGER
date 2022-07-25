# BIKER
"Bayesian Inference of the K600 Evasion Rate"

BIKER enables inference of the normalized riverine gas exchange velocity from measurements of river surface width and height/slope. It does so using Bayesian inference and a Hamiltonian Monte Carlo sampler to generate a posterior distribution for a gas exchange model. BIKER was developed in the context of the upcoming NASA/CNES/UKSA/CSA SWOT mission (https://swot.jpl.nasa.gov/mission/overview/).

## Installation
#### Dependencies
- R > 3.4.0
- methods
- Rcpp
- rstan
- rstantools
- dplyr
- reshape2
- rlang
- settings

#### To install
``` R
# First get devtools package
if (!require("devtools")) {
  install.packages("devtools")
  library("devtools")
}

#install BIKER
Sys.setenv(TAR = "/bin/tar") #if using a Linux machine (apparently)...
devtools::install_github("craigbrinkerhoff/BIKER", ref='main', force=TRUE)
```

## Getting started
See 'BIKER-manual.pdf' for all package information. Below are some notes to get you started on actually using the package.

To run BIKER, the following workflow is ideal. Note that BIKER follows a standard hydrology setup for the discretization of the river observations along a river reach. This means that matrix rows represent spatial steps along the reach while matrix columns represent timesteps for these observations.

#### Inputs
BIKER requires 3 inputs: <br>
- *Wobs*: a matrix of water surface widths <br>
- *Sobs*: a matrix of water surface slopes <br> 
- *dAobs*: a matrix of 'change in cross-sectional channel areas'. This must be approximated by the user. We generally assume a rectangular river channel and use the following R function to calculate dA from a timeseries of widths *w* and heights *h*.

``` R
#' Calculate partial cross-section area from width and height vectors (time series)
#' @param w vector of widths
#' @param h vector of heights(FROM MARK)
calcdA_vec <- function(w, h) {
  words <- order(w)
  warr <- w[words]
  harr <- h[words]
  delh <- c(0, diff(harr))
  delA <- cumsum(warr * delh)
  dA <- 1:length(w)
  dA[words] <- delA
  dA
}
```

#### Run BIKER
The following is the series of functions that need to be run (in this order) to use BIKER. Consult the help key in R to see examples. The below code will return the *k600* posterior mean and 95% confidence intervals (CIs) by default. The user can specify their CIs of choice in the biker_estimate function. Note that the 'meas_err' option should always be left to false as it is currently a work in progress and will produce erroneous *k600* estimates.<br>

``` R
reach_data <- biker_data(w=Wobs, s=Wobs, da=dAobs) #create an object of class bikerdata
reach_priors <- biker_priors(reach_data) #estimate prior hyperparameters for Bayesian inference using just river width and slope
k600_estimates <- biker_estimate(reach_data, reach_priors) #sample from joint posterior distribution to parameters
out <- biker_extract(k600_estimates) #extracts posterior means, sigmas, and CIs for parameters k600, n, and A0
```

#### Manually specifying prior hyperparameters
Priors can also be overwritten manually if you so desire. This is done as following, for example, for the A0_hat prior. Here, I set the prior to 5000 m2 and repeat for every spatial step in the river reach (because A0 is defined per spatial step):
``` R
reach_priors <- biker_priors(reach_data, logA0_hat = rep(5000, nrow(W_obs)))
```

If you want the full list of customizable prior hyperparameters, run the following:
``` R
prior_settings()
```

Finally, the 'cv2sigma()' function is included to if you want to specify hyperparameters by the coefficient of variation. This is mostly useful for the sigma hyperparameters.
