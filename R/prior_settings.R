#Default prior setup for algorithm





#' Options manager for BIKER default prior settings
#'
#' @param ... (Optional) named settings to query or set.
#' @param .__defaults See \code{?settings::option_manager}
#' @param .__reset See \code{?settings::option_manager}
#' @export
prior_settings <- settings::options_manager(
  paramnames = c("lowerbound_logk600", "upperbound_logk600", "lowerbound_A0",
                 "upperbound_A0", "lowerbound_logn", "upperbound_logn",
                 "sigma_post",
                 "logA0_hat", "logn_hat", "logk600_hat",
                 "logA0_sd", "logn_sd", "logk600_sd",
                 "River_Type", "k600_River_Type",
                 "Serr_sd", 'dAerr_sd'),

  # Bounds on parameters
  lowerbound_logk600 = log(0.001), #global-scope
  upperbound_logk600 = log(500), #global-scope

  lowerbound_A0 = rlang::quo(estimate_lowerboundA0(Wobs)),
  upperbound_A0 = rlang::quo(estimate_upperboundA0(Wobs)),
  lowerbound_logn = rlang::quo(estimate_lowerboundlogn(Wobs)),
  upperbound_logn = rlang::quo(estimate_upperboundlogn(Wobs)),

  # *Known* likelihood parameters
  sigma_post = 1.3, #mean sigma from 8,000 MC simulations of k600 model uncertainity

  # Hyperparameters via geoBAM & k600 prior
  logA0_hat = rlang::quo(estimate_logA0(Wobs)),
  logn_hat = rlang::quo(estimate_logn(Sobs, Wobs)),
  logk600_hat = rlang::quo(estimate_logk600(Wobs, Sobs)),

  #from geoBAM & k600 prior
  logA0_sd = rlang::quo(estimate_A0SD(Wobs)),
  logn_sd = rlang::quo(estimate_lognSD(Wobs)),
  logk600_sd = rlang::quo(estimate_logk600sd(Wobs, Sobs)),

  #Classified river type
  River_Type=rlang::quo(apply(Wobs, 1, classify_func)),
  k600_River_Type=rlang::quo(apply(Wobs, 2, classify_func_k600, Sobs=Sobs)),

  #SWOT Observation errors
  Serr_sd = 1.7e-5, #Durand et al. 2020 [unitless]: systematic SWOT error + layover error + random error
  dAerr_sd = rlang::quo(dA_sigma_func(Wobs)) #Durand et al. 2020 [m]
)
