#Functions to prepare SWOT observations for the inversion model


# Using advice from read-and-delete-me:
# "Be sure to add useDynLib(mypackage, .registration = TRUE) to the NAMESPACE file
# which you can do by placing the line   #' @useDynLib rstanarm, .registration = TRUE
# in one of your .R files
# see rstanarm's 'rstanarm-package.R' file"

#' Preprocess data for BIGEER estimation
#'
#' Produces a bigeedata object that can be passed to bigee_estimate function
#'
#' @useDynLib BIGEER, .registration = TRUE
#' @param w Matrix (or data frame) of widths: time as columns, space as rows
#' @param s Matrix of slopes: time as columns, space as rows
#' @param dA Matrix of area above base area: time as columns, space as rows
#' @param max_xs Maximum number of cross-sections to allow in data. Used to reduce
#'   sampling time. Defaults to 30.
#' @param seed RNG seed to use for sampling cross-sections, if nx > max_xs.
#' @export

bigee_data <- function(w,
                     s,
                     dA,
                     max_xs = 30L,
                     seed = NULL) {

  manning_ready <- !is.null(s) && !is.null(dA)
  if (!manning_ready) {
    s <- dA <- matrix(1, nrow = nrow(w), ncol = ncol(w))
  }


  datalist <- list(Wobs = w,
                Sobs = s,
                dAobs = dA)

  datalist <- bigee_check_args(datalist)
  datalist <- bigee_check_nas(datalist)

  nx <- nrow(datalist$Wobs)
  nt <- ncol(datalist$Wobs)

  out <- structure(c(datalist,
                        nx = nx,
                        nt = nt),
                   manning_ready = manning_ready,
                   class = c("bigeedata"))

  if (nx > max_xs)
    out <- sample_xs(out, n = max_xs, seed = seed)

  out
}

#' Performs the following checks:
#' - types:
#'     - everything else matrix
#' - dimensions:
#'     - all matrices have same dims
#'
#' @param datalist A list of BIGEE data inputs
bigee_check_args <- function(datalist) {
  matlist <- datalist[names(datalist)]

  if (!all(vapply(matlist, is, logical(1), "matrix")))
    stop("All data must be a supplied as a matrix.\n")


  # Check dims
  nr <- nrow(matlist[[1]])
  nc <- ncol(matlist[[1]])
  if (!(all(vapply(matlist, nrow, 0L) == nr) &&
        all(vapply(matlist, ncol, 0L) == nc)))
    stop("All data must have same dimensions.\n")

  out <- c(matlist)

  out
}

#' Add missing-data inputs to data list
#'
#' Binary matrices indicating where data are/aren't missing are
#' added to the data list. This is required in order to run
#' ragged-array data structures in the stanfile.
#'
#' Previously this function omitted any times with missing data,
#' but now that ragged arrays are accommodated in the stanfile the
#' operations are entirely different.
#'
#' @param datalist a list of BIGEER inputs
#' @importFrom stats median
bigee_check_nas <- function(datalist) {

  mats <- vapply(datalist, is.matrix, logical(1))
  nonas <- lapply(datalist[mats], function(x) !is.na(x))

  # Manning has-data matrix (only nonzero if all Manning obs present)
  if (identical(setdiff(c("Wobs", "Sobs", "dAobs"),
                        names(datalist[mats])),
                character(0))) {
    hasdat_s <- (!is.na(datalist[["Sobs"]])) * 1
    hasdat_a <- (!is.na(datalist[["dAobs"]])) * 1
    hasdat_w <- (!is.na(datalist[["Wobs"]])) * 1

    hasdat_man <- hasdat_w * hasdat_s * hasdat_a

    # Replace NA's with zeros so Stan will accept the data
    datalist[["Sobs"]][!hasdat_man] <- 0
    datalist[["dAobs"]][!hasdat_man] <- 0
    datalist[["Wobs"]][!hasdat_man] <- 0
  } else {
    stop('Not all Mannings observables are present!')
  }

  #check that there are at least 3 non-zero values in each row and/or column
  w_test_col <- apply(datalist[["Wobs"]],2, function(x){sum(x>0, na.rm=T)})  #make sure at least 3 observations for width sd calculation for priors
  w_test_row <- apply(datalist[["Wobs"]],1, function(x){sum(x>0, na.rm=T)})  #make sure at least 3 observations for width sd calculation for priors
  if(any(w_test_col<2 | w_test_row < 2)){
    stop('Need at least 3 non-NA values in each row and column to calculate some priors!')}

  #dA error checks
  if (!is.null(datalist[["dAobs"]])) {
    dA_shift <- apply(datalist[["dAobs"]], 1, function(x) median(x) - min(x))
  } else {
    dA_shift <- rep(0, nrow(datalist[["Wobs"]]))
  }

  newbits <- list(
    hasdat_man = hasdat_man,
    ntot_man = sum(hasdat_man),
    dA_shift = dA_shift
  )

  out <- c(datalist, newbits)
  out
}

#' Establish prior hyperparameters for BIGEER estimation
#'
#' Produces a bigeepriors object that can be passed to bigee_estimate function
#'
#' @useDynLib BIGEER, .registration = TRUE
#' @param bigeedata An object of class bigeedata, as returned by \code{bigee_data}
#' @param ... Optional manually set parameters. Unquoted expressions are allowed,
#'   e.g. \code{logk600_sd = cv2sigma(0.8)}. Additionally, any variables present in
#'   \code{bigeedata} may be referenced, e.g. \code{lowerbound_logk600 = log(mean(Wobs)) + log(5)}
#' @export

bigee_priors <- function(bigeedata,
                       ...) {
  force(bigeedata)
  paramset <- prior_settings("paramnames")

  myparams0 <- rlang::quos(..., .named = TRUE)
  myparams <- do.call(settings::clone_and_merge,
                      args = c(list(options = prior_settings), myparams0))

  quoparams <- myparams()[-1] # first one is parameter set
  params <- lapply(quoparams, rlang::eval_tidy, data = bigeedata)

  if (!length(params[["logk600_sd"]]) == bigeedata$nt)
    params$logk600_sd <- rep(params$logk600_sd, length.out = bigeedata$nt)

  if (!identical(dim(params[["sigma_post"]]),
                 as.integer(c(bigeedata$nx, bigeedata$nt)))) {
    params$sigma_post <- matrix(rep(params$sigma_post,
                                   length.out = bigeedata$nt * bigeedata$nx),
                               nrow = bigeedata$nx)
  }

  #just priors the user wants to see
  sigma_paramset <- c('sigma_post')
  sigma_paramset <- params[sigma_paramset]

  #total priors needed to run geoBAM
  bigee_paramset <- c("lowerbound_logk600", "upperbound_logk600", "lowerbound_A0",
                      "upperbound_A0", "lowerbound_logn", "upperbound_logn",
                      "logA0_hat", "logn_hat","logk600_hat", "logA0_sd", "logn_sd", "logk600_sd")
  bigeeparams <- params[bigee_paramset]

  riverType <- params[c("River_Type", "k600_River_Type")]

  out <- list( 'river_types'=riverType, 'river_type_priors'=bigeeparams, 'sigma_model'=sigma_paramset)
  out <- structure(out,
                   class = c("bigeepriors"))
  out
}

compose_bigee_inputs <- function(bigeedata, priors = bigee_priors(bigeedata)) {

  inps <- c(bigeedata, priors)

  out <- inps
  out

}


#' Take a random sample of a bigeedata object's cross-sections.
#'
#' @param bigeedata a bigeedata object, as returned by \code{bigee_data()}
#' @param n Number of cross-sections to
#' @param seed option RNG seed, for reproducibility.
#' @importFrom methods is
#' @export
sample_xs <- function(bigeedata, n, seed = NULL) {

  stopifnot(is(bigeedata, "bigeedata"))

  if (n >= bigeedata$nx)
    return(bigeedata)

  if (!is.null(seed))
    set.seed(seed)
  keepxs <- sort(sample(1:bigeedata$nx, size = n, replace = FALSE))

  bigeedata$nx <- n
  bigeedata$Wobs <- bigeedata$Wobs[keepxs, ]

  if (!is.null(bigeedata$Sobs)) {
    bigeedata$Sobs <- bigeedata$Sobs[keepxs, ]
    bigeedata$dAobs <- bigeedata$dAobs[keepxs, ]
  }

  bigeedata
}



#' Calculate lognormal moments based on truncated normal parameters
#'
#' Used to put measurement errors into original log-normal parameterization.
#'
#' @param obs A numeric vector of observations
#' @param err_sigma Standard deviation of measurement error
#' @param a zero-reference for method of moments.
#' @importFrom stats dnorm pnorm
ln_moms <- function(obs, err_sigma, a = 0) {
  alpha <- (a - obs) / err_sigma
  Z <- 1 - pnorm(alpha)

  mean <- obs + (dnorm(alpha)) / Z * err_sigma
  sdquan <- 1 + (alpha * dnorm(alpha)) / Z -
    (dnorm(alpha) / Z)^2
  sd <- err_sigma * sqrt(sdquan)

  out <- list(mean = mean, sd = sd)
  out
}

#' Calculate lognormal sigma parameter based on truncated normal parameters
#'
#' Used to put measurement errors into original log-normal parameterization.
#'
#' @param obs A numeric vector of observations
#' @param err_sigma Standard deviation of measurement error
#' @param a zero-reference for method of moments.
ln_sigsq <- function(obs, err_sigma, a = 0) {
  moms <- ln_moms(obs = obs, err_sigma = err_sigma, a = a)
  mn <- unname(moms[["mean"]])
  sd <- unname(moms[["sd"]])
  mu <- 2 * log(mn) - 0.5 * log(sd^2 + mn^2)
  sigsq <- 2 * log(mn) - 2 * mu

  sigsq
}
