
#' Options manager for BIKER default prior settings for k600 model
#'
#' @param ... (Optional) named settings to query or set.
#' @param .__defaults See \code{?settings::option_manager}
#' @param .__reset See \code{?settings::option_manager}
#' @export
prior_settings <- settings::options_manager(
  paramnames = c("lowerbound_logk", "upperbound_logk", "lowerbound_A0", "upperbound_A0", "lowerbound_logn", "upperbound_logn",
                 "sigma_post",
                 "logA0_hat", "logk_hat", "logn_hat",
                 "logA0_sd", "logk_sd", "logn_sd",
                 "River_Type",
                 "Serr_sd", 'dAerr_sd'),

  # Bounds on parameters
  lowerbound_logk = log(0.001), #global-scope
  upperbound_logk = log(500), #global-scope

  lowerbound_A0 = rlang::quo(estimate_lowerboundA0(Wobs)),
  upperbound_A0 = rlang::quo(estimate_upperboundA0(Wobs)),

  lowerbound_logn = log(0.01),
  upperbound_logn = log(0.05),

  # *Known* likelihood parameters
  sigma_post = 1.12, #standard error of k~ustar model parameter

  # Hyperparameters via geoBAM & k prior
  logA0_hat = rlang::quo(estimate_logA0(Wobs)),
  logk_hat = rlang::quo(estimate_logk(Sobs)),
  logn_hat = rlang::quo(estimate_logn(Wobs, Sobs)),

  #from geoBAM & k prior
  logA0_sd = rlang::quo(estimate_A0SD(Wobs)),
  logk_sd = rlang::quo(estimate_logksd(Sobs)),
  logn_sd = rlang::quo(estimate_lognSD(Wobs)),

  #Classified river type
  River_Type=rlang::quo(apply(Wobs, 1, classify_func)),

  #SWOT Observation errors following Hagemann etal 2017
  Serr_sd = 1e-5,
  dAerr_sd = 10,
  Werr_sd = 10
)

