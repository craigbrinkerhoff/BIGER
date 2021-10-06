#Functions to prepare SWOT observations for the inversion model


# Using advice from read-and-delete-me:
# "Be sure to add useDynLib(mypackage, .registration = TRUE) to the NAMESPACE file
# which you can do by placing the line   #' @useDynLib rstanarm, .registration = TRUE
# in one of your .R files
# see rstanarm's 'rstanarm-package.R' file"

#' Preprocess data for BIKER estimation
#'
#' Produces a bikerdata object that can be passed to biker_estimate function
#'
#' @useDynLib BIKER, .registration = TRUE
#' @param w Matrix (or data frame) of widths: time as columns, space as rows
#' @param s Matrix of slopes: time as columns, space as rows
#' @param dA Matrix of area above base area: time as columns, space as rows
#' @param priorQ Mean annual flow prior for Q
#' @param max_xs Maximum number of cross-sections to allow in data. Used to reduce
#'   sampling time. Defaults to 30.
#' @param seed RNG seed to use for sampling cross-sections, if nx > max_xs.
#' @export

biker_data <- function(w,
                     s,
                     dA,
                     priorQ,
                     max_xs = 30L,
                     seed = NULL) {

  manning_ready <- !is.null(s) && !is.null(dA)
  if (!manning_ready) {
    s <- dA <- matrix(1, nrow = nrow(w), ncol = ncol(w))
  }


  datalist <- list(Wobs = w,
                Sobs = s,
                dAobs = dA,
                priorQ=priorQ)

  datalist <- biker_check_args(datalist)
  datalist <- biker_check_nas(datalist)

  nx <- nrow(datalist$Wobs)
  nt <- ncol(datalist$Wobs)

  out <- structure(c(datalist,
                        nx = nx,
                        nt = nt),
                   manning_ready = manning_ready,
                   class = c("bikerdata"))

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
#' @param datalist A list of biker data inputs
biker_check_args <- function(datalist) {
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
#' @param datalist a list of BIKER inputs
#' @importFrom stats median
biker_check_nas <- function(datalist) {

  mats <- vapply(datalist, is.matrix, logical(1))
  nonas <- lapply(datalist[mats], function(x) !is.na(x))

  # Manning has-data matrix (only nonzero if all Manning obs present)
  if (identical(setdiff(c("Wobs", "Sobs", "dAobs"),
                        names(datalist[mats])),
                character(0))) {
    hasdat_s <- (!is.na(datalist[["Sobs"]])) * 1
    hasdat_a <- (!is.na(datalist[["dAobs"]])) * 1
    hasdat_w <- (!is.na(datalist[["Wobs"]])) * 1

    hasdat <- hasdat_w * hasdat_s * hasdat_a

    # Replace NA's with zeros so Stan will accept the data
    datalist[["Sobs"]][!hasdat] <- 0
    datalist[["dAobs"]][!hasdat] <- 0
    datalist[["Wobs"]][!hasdat] <- 0
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
    hasdat = hasdat,
    ntot_man = sum(hasdat),
    dA_shift = dA_shift
  )

  out <- c(datalist, newbits)
  out
}

#' Establish prior hyperparameters for BIKER estimation
#'
#' Produces a bikerpriors object that can be passed to biker_estimate function
#'
#' @useDynLib BIKER, .registration = TRUE
#' @param bikerdata An object of class bikerdata, as returned by \code{biker_data}
#' @param ... Optional manually set parameters. Unquoted expressions are allowed,
#'   e.g. \code{logk_sd = cv2sigma(0.8)}. Additionally, any variables present in
#'   \code{bikerdata} may be referenced, e.g. \code{lowerbound_logk = log(mean(Wobs)) + log(5)}
#' @export
biker_priors <- function(bikerdata,
                        ...) {
  force(bikerdata)
  paramset <- prior_settings("paramnames")

  myparams0 <- rlang::quos(..., .named = TRUE)
  myparams <- do.call(settings::clone_and_merge,
                      args = c(list(options = prior_settings), myparams0))

  quoparams <- myparams()[-1] # first one is parameter set
  params <- lapply(quoparams, rlang::eval_tidy, data = bikerdata)

  if (!length(params[["logk_sd"]]) == bikerdata$nt)
    params$logk6_sd <- rep(params$logk_sd, length.out = bikerdata$nt)

  if (!identical(dim(params[["sigma_post"]]),
                 as.integer(c(bikerdata$nx, bikerdata$nt)))) {
    params$sigma_post <- matrix(rep(params$sigma_post,
                                   length.out = bikerdata$nt * bikerdata$nx),
                               nrow = bikerdata$nx)
  }

  #just priors the user wants to see
  sigma_paramset <- c('sigma_post', 'Serr_sd', 'dAerr_sd', 'Werr_sd')
  sigma_paramset <- params[sigma_paramset]

  #total priors needed to run geoBAM
  biker_paramset <- c("lowerbound_logk", "upperbound_logk", "lowerbound_A0","upperbound_A0", "lowerbound_logn", "upperbound_logn",
                      "logA0_hat","logk_hat","logn_hat",
                      "logA0_sd", "logk_sd", "logn_sd")
  bikerparams <- params[biker_paramset]

  riverType <- params[c("River_Type")]

  out <- list( 'river_types'=riverType, 'river_type_priors'=bikerparams, 'sigma_model'=sigma_paramset)
  out <- structure(out,
                   class = c("bikerpriors"))
  out
}

compose_biker_inputs <- function(bikerdata, priors = biker_priors(bikerdata)) {

  inps <- c(bikerdata, priors)

  out <- inps
  out

}


#' Take a random sample of a bikerdata object's cross-sections.
#'
#' @param bikerdata a bikerdata object, as returned by \code{biker_data()}
#' @param n Number of cross-sections to
#' @param seed option RNG seed, for reproducibility.
#' @importFrom methods is
#' @export
sample_xs <- function(bikerdata, n, seed = NULL) {

  stopifnot(is(bikerdata, "bikerdata"))

  if (n >= bikerdata$nx)
    return(bikerdata)

  if (!is.null(seed))
    set.seed(seed)
  keepxs <- sort(sample(1:bikerdata$nx, size = n, replace = FALSE))

  bikerdata$nx <- n
  bikerdata$Wobs <- bikerdata$Wobs[keepxs, ]

  if (!is.null(bikerdata$Sobs)) {
    bikerdata$Sobs <- bikerdata$Sobs[keepxs, ]
    bikerdata$dAobs <- bikerdata$dAobs[keepxs, ]
  }

  bikerdata
}
