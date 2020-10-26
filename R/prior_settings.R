#' Options manager for BIGEER default prior settings
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
                 "River_Type", "k600_River_Type"),

  # Bounds on parameters
  lowerbound_logk600 = rlang::quo(estimate_lowerboundlogk600(Wobs)),
  upperbound_logk600 = rlang::quo(estimate_upperboundlogk600(Wobs)),

  lowerbound_A0 = rlang::quo(estimate_lowerboundA0(Wobs)),
  upperbound_A0 = rlang::quo(estimate_upperboundA0(Wobs)),
  lowerbound_logn = rlang::quo(estimate_lowerboundlogn(Wobs)),
  upperbound_logn = rlang::quo(estimate_upperboundlogn(Wobs)),

  # *Known* likelihood parameters
  sigma_post = 0.30, #obtained from monte carlo sampling of k600 model

  # Hyperparameters via geoBAM & k600 prior
  logA0_hat = rlang::quo(estimate_logA0(Wobs)),
  logn_hat = rlang::quo(estimate_logn(Sobs, Wobs)),
  logk600_hat = rlang::quo(estimate_logk600(Wobs)),

  #from geoBAM & k600 prior
  logA0_sd = rlang::quo(estimate_A0SD(Wobs)),
  logn_sd = rlang::quo(estimate_lognSD(Wobs)),
  logk600_sd = rlang::quo(estimate_logk600sd(Wobs)),

  #Classified river type
  River_Type=rlang::quo(apply(Wobs, 1, classify_func)),
  k600_River_Type=rlang::quo(apply(Wobs, 2, classify_func_k600))
)
