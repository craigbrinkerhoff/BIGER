# BIKER
"Bayesian Inference/Inversion of the K600 Evasion Rate"

BIKER enables estimation of the normalized riverine gas exchange velocity from measurements of river surface width and height/slope. It does so using Bayesian inference and a Hamiltonian Monte Carlo sampler to generate a posterior distribution for a gas exchange model. BIKER was developed in the context of the NASA/CNES/UKSA/CSA SWOT mission (https://swot.jpl.nasa.gov/mission/overview/). The academic manuscript associated with this model development and validation is in preperation. The code for developing and validing BIKER (as well as writing the manuscript in RMarkdown) is available at https://github.com/craigbrinkerhoff/RSK600.

## Installation
```
Sys.setenv(TAR = "/bin/tar") #if using a Linux machine...
devtools::install_github("craigbrinkerhoff/BIKER", ref='main', force=TRUE)
```

## Getting started
To run BIKER, the following workflow is ideal. Note that BIKER follows a standard hydrology setup for the discretization of the river observations along a mass-conserved river reach. This means that matrix rows represent spatial steps along the reach while matrix columns represent timesteps for these observations.

#### Inputs
BIKER requires 3 inputs: <br>
- $Wobs$: a matrix of water surface widths <br>
- $Sobs$: a matrix of water surface slopes <br> 
- $dAobs$: a matrix of 'change in cross-sectional channel areas'. This must be approximated by the user. We generally assume a rectangular river channel so that this can be calculated as the following: $Wobs*\delta Hobs$. $Hobs$ is the matrix of water-surface elevations.

#### Run BIKER
The following is the series of fuctions that need to be run (in this order) to use BIKER. Consult the help key in R to see examples. The below code will return the $k_{600}$ posterior mean and 95% confidence intervals (CIs) by default. The user can specify their CIs of choice in the biker_estimate function. Note that the 'meas_err' option should always be left to false as it is currently a work in progress and will produce erronous $k_{600}$ estimates.<br>

```
reach_data <- biker_data(w=Wobs, s=Wobs, da=dAobs) #collect data into object the algorithm can read
reach_priors <- biker_priors(reach_data) #estimate prior hyperparameters for Bayesian inference using just river width and slope
k600_estimates <- biker_estimate(reach_data, reach_priors) #sample from joint posterior distribution to obtain estimates of k600
```

#### Manually specifying prior hyperparameters
Priors can also be overwritten manually if you so desire. This is done as following, for example, for the A0_hat prior. Here, I set the prior to 5000 m2 and repeat for every spatial step in the river reach (because A0 is defined per spatial step):
```
reach_priors <- biker_priors(reach_data, logA0_hat = rep(5000, nrow(W_obs)))
```

If you want the full list of customizable prior hyperparameters, run the following:
```
prior_settings()
```

Finally, the 'cv2sigma()' function is included to if you want to specify hyperparameters by the cofficient of variation. This is mostly useful for the sigma hyperparameters.
