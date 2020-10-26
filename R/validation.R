# validation functions

## Model performance metrics --------------------

# Sum of squares function
sumsq <- function(x) sum(x^2)

#' Create a data.frame for BAM validation
#'
#' @param fit A stanfit object, as returned from \code{bigee_estimate()}
#' @param qobs a vector of observed flow.
#'
#' @importFrom rlang .data
#' @export
bigee_valdata <- function(fit, k600obs) {
  stopifnot(is(fit, "stanfit"))
  stopifnot(is.numeric(k600obs))
  k600pred <- bigee_k600pred(fit = fit, chain = "all") %>%
    dplyr::transmute(.data$time, k600pred = mean)
  stopifnot(length(k600obs) == nrow(k600pred))
  out <- cbind(k600pred, k600obs = k600obs)
  out
}

#' Calculate validation metrics and plots
#'
#' @param fit A stanfit object, as returned from \code{bigee_estimate()}
#' @param k600obs a vector of observed flow.
#' @param stats Which stats to include in the summary?
#'
#' @export
bigee_validate <- function(fit, k600obs, stats = c("RRMSE", "MRR", "SDRR",
                                              "NRMSE", "rBIAS",
                                              "CoV", 'r2')) {
  stats <- match.arg(stats, several.ok = TRUE)
  valdata <- bigee_valdata(fit = fit, k600obs = k600obs)
  pred <- valdata$k600pred
  obs <- valdata$k600obs

  statvals <- vapply(stats, do.call, numeric(1),
                     args = list(pred = pred, meas = obs))

  out <- structure(list(valdata = valdata,
                        stats = statvals),
                   class = c("bigeeval"))
}


#' Relative root-mean-square error
#'
#' @param pred vector of predictions
#' @param meas vector of measurements
#' @export
RRMSE <- function(pred, meas)
  sqrt(mean((pred - meas)^2 / meas^2))

#' Mean relativ residual
#'
#' @param pred vector of predictions
#' @param meas vector of measurements
#' @export
MRR <- function(pred, meas)
  mean((meas - pred) / meas)

#' Standard deviation of relative residual
#'
#' @param pred vector of predictions
#' @param meas vector of measurements
#' @export
SDRR <- function(pred, meas)
  sd((meas - pred) / meas)


#' Normalized root-mean-square error
#'
#' @param pred vector of predictions
#' @param meas vector of measurements
#' @export
NRMSE <- function(pred, meas)
  sqrt(mean((meas - pred)^2)) / mean(meas)

#' Relative bias
#'
#' @param pred vector of predictions
#' @param meas vector of measurements
#' @export
rBIAS <- function(pred, meas)
  mean(pred - meas) / mean(meas)

#' Coefficient of variation
#'
#' @param pred vector of predictions
#' @param meas vector of measurements
#' @importFrom stats sd
#' @export
#'
CoV <- function(pred, meas)
  sd(pred - meas) / mean(meas)

#' Coefficient of determination
#'
#' @param pred vector of predictions
#' @param meas vector of measurements
#' @importFrom stats sd
#' @export
r2 <- function(pred, meas){
  1 - (sumsq(meas-pred)/sumsq(meas-mean(meas)))
}

