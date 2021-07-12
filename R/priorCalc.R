#Functions for default prior setup for algorithm (n and A0 functions from Brinkerhoff etal 2020 --> geoBAM-Expert algorithm)



#k prior calculations--------------------------------------------------------------------------------------------------
#' Estimate k_hat using biker data for k600 model
#'
#' @param Sobs Observed S,as a space-down, time-across matrix
#' @export
estimate_logk <- function(Sobs){
  Sobs[Sobs <=0] <- NA

  colSobs <- colMeans(log(Sobs), na.rm=T)
 # khat <- 5.0941 + 0.6417*colSobs #ifelse(colSobs < -4.634, 3.22 + 0.347*colSobs, 6.85 + 1.13*colSobs)
  khat <- log(1144.6*colMeans(Sobs, na.rm=T)) #r2 0.52 using ulseth data where Rh=H
}

#' Estimate k sd prior using biker data for k600 model
#'
#' @param Sobs Observed S,as a space-down, time-across matrix
#' @export
estimate_logksd <- function(Sobs){
  Sobs[Sobs <=0] <- NA

  ksd <- rep(3.841, ncol(Sobs)) #standard model error of k600=beta*slope model fit to Ulseth data
}

# Prior calculation using geoBAM-Expert classification framework------------------------------------------------------------------
#class 17 are 'big' riverrs

#' Estimate base cross-sectional area using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix
#' @export
estimate_logA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #Expert classification
  temp <- c(3.723280881,
            4.497073804,
            4.753590191,
            4.990432587,
            4.911328517,
            5.350615603,
            5.422742508,
            5.523458921,
            5.774532256,
            6.446836611,
            6.527953643,
            6.873163834,
            7.102499356,
            8.007965013,
            8.937204637,
            4.432006567)

  class <- apply(Wobs, 1, classify_func)
  logA0_hat <- ifelse(class != 17, temp[class], -0.2918 + 1.6930 * lwbar - 0.1887 * lwsd) #repeat for each sptial unit
  #global r2: 0.907
}

#' Estimate base cross-sectional area lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification
  temp <- c(0.732367894,
            1.508511994,
            0.91027266,
            2.517696473,
            1.199964783,
            2.681021529,
            2.148850993,
            1.545432582,
            2.415913778,
            3.106826321,
            3.874321138,
            2.694627181,
            3.696351469,
            3.593194204,
            4.043928076,
            0.262364264)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_A0 <- ifelse(class != 17, exp(temp[class]), exp(4.540631665))
  lowerbound_A0 <- min(lowerbound_A0, na.rm = TRUE)
}

#' Estimate base cross-sectional area upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification
  temp <- c(7.640123173,
            7.355641103,
            8.997147152,
            9.164296433,
            8.554488976,
            9.417354541,
            7.677863501,
            8.144679183,
            7.863266724,
            8.793308627,
            8.776475789,
            9.014325488,
            8.78186249,
            9.61580548,
            11.6483301,
            11.55214618)

  class <- apply(Wobs, 1, classify_func)
  upperbound_A0 <- ifelse(class != 17, exp(temp[class]), 1000000)
  upperbound <- max(upperbound_A0, na.rm = TRUE)
}

#' Estimate base cross-sectional area SD using bam dat
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_A0SD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(1.446820524,
            1.350055707,
            1.363875688,
            1.260787706,
            1.306619211,
            1.315826017,
            1.190542029,
            1.337220271,
            1.21047899,
            1.359096608,
            1.245863462,
            1.287297345,
            1.08535437,
            1.154319081,
            1.5575699,
            2.272342301)

  class <- apply(Wobs, 1, classify_func)
  logA0_sd <- ifelse(class != 17, temp[class], 0.58987527)
}

