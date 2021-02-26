# BIKER
"Bayesian Inference/Inversion of the K600 Evasion Rate"

BIKER uses Bayesian inference to simultaneously estimate river channel size, roughness, and the normalized gas exchange velocity from some mass-conserved river reach. BIKER was developed specifically to ingest observations from the NASA/CNES/UKSA/CSA Surface Water and Ocean Topography (SWOT) mission, but can theoreitcally be used on any dataset of river width and water surface elevation measurements made along a set of mass-conserved river cross-sections.

The academic manuscript associated with BIKER and its validation is in preperation. The code for developing and validing BIKER and writing the manuscript is available at https://github.com/craigbrinkerhoff/RSK600.


## Installation
```
Sys.setenv(TAR = "/bin/tar") #if using a Linux machine...
devtools::install_github("craigbrinkerhoff/BIKER", ref='main', force=TRUE)
```

## Getting started
BIKER requires three inputs to the algorithm for some mass-conserved river reach: a matrix of river width measurements $W_{obs}$, a matrix of water surface elevation measurements $S_{obs}$, and a matrix of 'change in cross-sectional channel area' observations $dA_{obs}$. All three BIKER inputs must follow the standard hydrology discretization setup where spatial steps are row-wise and temporal steps are column-wise. $dA_{obs}$ can be estimated assuming a rectangular river channel using the following expression: $dA_{obs}\approx W_{obs}* \delta H_{obs}$ for water surface elevation $H_{obs}$.

BIKER requires the user to run three seperate functions in order to estimate $k_{600}$. First, you must gather all necessary data into a BIKER object that the sampling model can ingest. Then, you need to set values for your prior hyperparameters, and then you can "turn the Bayesian crank" to obtain $k_{600}$ estimates. See the code below for doing this.

```
#set up biker object for some mass conserved river reach
data <- biker_data(w=W_obs, s=S_obs, dA=dA_obs)

#Use internal functions to automatically assign hyperparameters. These can be overwritten manually too
prior_hyperparameters <- biker_priors(data)

#run BIKER
k600_estimate <- biker_estimate(biker_data, prior_hyperparameters)
```
This will return the $k_{600}$ posterior mean and 95% confidence intervals (CIs) by default. The user can specify their CIs of choice in the biker_estimate function. Note that the 'meas_err' option should always be left to false as it is currently a work in progress and will produce erronous $k_{600}$ estimates.

Also note that to manually override the hyperparameter assignments, one simply has to specify the specific prior they want to override. All geomorphic priors are within the 'river_type_priors' sublist. For example, lets use the upperbound on Manning's n:
```
prior_hyperparameters <- biker_priors(data)
prior_hyperparameters$river_type_priors$upperbound_logn = -2

#one can also change the upscaling model uncertainity as follows
prior_hyperparameters$sigma_model$sigma_post = 0.25

#for the full list of hyperparameters that can be changed, run this:
prior_settings()
```

## Troubleshooting
Use the help key in RStudio to get function options and email cbrinkerhoff[at]umass[dot]edu for further help!
