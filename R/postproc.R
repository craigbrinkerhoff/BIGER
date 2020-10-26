# Functions to aid in post-processing BIGGE results

#' k600 posterior mean and Bayesian credible interval.
#' 
#' @param fit A stanfit object, as returned from \code{bigge_estimate()}
#' @param chain Either an integer specifying which chain(s) to extract statistics from,
#'   or "all" (the default), in which case they are extracted from all chains.
#' @param conf.level A numeric value on (0, 1) specifying the size of the Bayesian 
#'   credible interval. Default is 0.95.
#' @importFrom stats quantile
#' @importFrom rlang .data
#' @export 

bigge_k600pred <- function(fit, chain = "all", conf.level = 0.95) {
  
  k600post <- rstan::extract(fit, "logk600", permuted = FALSE) %>% 
    reshape2::melt()
  
  if (conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be on the interval (0,1).\n")
  
  alpha <- 1 - conf.level
  
  nchains <- fit@sim$chains
  if (chain == "all") 
    chain <- 1:nchains
  stopifnot(is.numeric(chain))
  
  k600stats <- k600post %>% 
    dplyr::mutate(chains = gsub("^chain:", "", .data$chains)) %>% 
    dplyr::filter(.data$chains %in% chain) %>% 
    dplyr::mutate(value = exp(.data$value)) %>% 
    dplyr::group_by(.data$parameters) %>%
    dplyr::summarize(mean = mean(.data$value),
              conf.low = quantile(.data$value, alpha / 2),
              conf.high = quantile(.data$value, 1 - (alpha / 2))) %>% 
    dplyr::rename(time = .data$parameters) %>% 
    dplyr::mutate(time = gsub("^logk600\\[", "", .data$time),
           time = gsub("\\]$", "", .data$time),
           time = as.numeric(.data$time)) %>% 
    dplyr::arrange(.data$time)
  
  k600stats
}
