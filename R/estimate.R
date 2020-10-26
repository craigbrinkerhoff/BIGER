
#' Estimate BIGGE
#'
#' Fits a BIGGE model of one of several variants using Hamiltonian Monte Carlo.
#'
#' @param biggedata A bamdata object, as produced by \code{bigge_data()}
#' @param biggepriors A bampriors object. If none is supplied, defaults are used
#'   from calling \code{bigge_priors(biggedata)} (with no other arguments).
#' @param cores Number of processing cores for running chains in parallel.
#'   See \code{?rstan::sampling}. Defaults to \code{parallel::detectCores()}.
#' @param chains A positive integer specifying the number of Markov chains.
#'   The default is 3.
#' @param iter Number of iterations per chain (including warmup). Defaults to 1000.
#' @param stanmodel A \code{stanmodel} object to use instead of one of the default
#'   models.
#' @param pars (passed to \code{rstan::sampling()}) A vector of character strings specifying
#'   parameters of interest to be returned in the stanfit object. If not specified,
#'   a default parameter set is returned.
#' @param include (passed to \code{rstan::sampling()}) Defaults to FALSE, which
#'   excludes parameters specified in \code{pars} from the returned model.
#' @param ... Other arguments passed to rstan::sampling() for customizing the
#'   Monte Carlo sampler
#' @import rstan
#' @export

bigge_estimate <- function(biggedata,
                         biggepriors = NULL,
                         cores = getOption("mc.cores", default = parallel::detectCores()),
                         chains = 3L,
                         iter = 1000L,
                         stanmodel = NULL,
                         pars = NULL,
                         include = FALSE,
                         ...) {
  stopifnot(is(biggedata, "biggedata"))
  if (is.null(biggepriors))
    biggepriors <- bigge_priors(biggedata)
  stopifnot(is(biggepriors, "biggepriors"))

  #Reforamte priors to a single list for stan
  biggepriors <- c(biggepriors[[2]], biggepriors[[3]])

  biggeinputs <- compose_bigge_inputs(biggedata, biggepriors)
  biggeinputs$inc_m <- 1

  if (!is.null(stanmodel)) {
    stopifnot(inherits(stanmodel, "stanmodel"))
    stanfit <- stanmodel
  } else {
    stanfit <- stanmodels[["master"]]
  }

  if (is.null(pars)) {
    pars <- c("man_rhs", "logWSpart",
              "logk600tn", "logk600nbar",
              "Sact", "Wact", "dAact")
  }

  out <- sampling(stanfit, data = biggeinputs,
                  cores = cores, chains = chains, iter = iter,
                  pars = pars, include = include,
                  ...)

  out
}


