#Functions for default prior setup for algorithm (n and A0 functions from Brinkerhoff etal 2020 --> geoBAM-Expert algorithm)



#k600 prior calculations--------------------------------------------------------------------------------------------------
#priors are based on width classes from Ulseth et al. (2019) field dataset of k600


#' Estimate k600_hat using bigee data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix
#' @export
estimate_logk600 <- function(Wobs, Sobs){
  Wobs[Wobs <= 0] <- NA
  Sobs[Sobs <=0] <- NA
  temp <- log(c(5.666775,
                3.024020,
                2.034140,
                2.238760,
                79.464785))

  class <- apply(Wobs, 2, classify_func_k600, Sobs=Sobs)
  khat <- temp[class]
}

#' Estimate k600 sd prior using bigee data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix
#' @export
estimate_logk600sd <- function(Wobs, Sobs){
  Wobs[Wobs <= 0] <- NA
  Sobs[Sobs <=0] <- NA
  temp <- log(c(20.274259,
                16.988323,
                4.536993,
                3.044974,
                653.966672))

  class <- apply(Wobs, 2, classify_func_k600, Sobs=Sobs)
  ksd <- temp[class]
}


# Prior calculation using expert classification framework------------------------------------------------------------------
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

#'Estimate manning's n using bam data
#'
#' @param Sobs Observed S, as a space-down, time-across matrix
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_logn <- function(Sobs, Wobs) {
  Sobs[Sobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lsbar <- apply(log(Sobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(-2.956882373,
            -3.218074003,
            -3.073438472,
            -3.316850463,
            -3.318215717,
            -3.15335243,
            -3.456243594,
            -3.545592423,
            -3.240085716,
            -3.401877538,
            -3.269861085,
            -3.372948626,
            -3.404202242,
            -3.274729636,
            -3.405497386,
            -3.318138928)

  class <- apply(Wobs, 1, classify_func)
  logn_hat <- ifelse(class != 17, temp[class], -0.1636 + 0.4077 * lsbar) #repeat for each sptial unit
  #Global r2: 0.631
}

#'Estimate manning's n lowerbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lowerboundlogn <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(-5.678333858,
            -5.467593416,
            -5.67856698,
            -5.768715685,
            -5.682940688,
            -5.841321589,
            -5.56434293,
            -5.957755533,
            -5.646550723,
            -6.064509006,
            -6.029056118,
            -6.405084634,
            -6.678294272,
            -6.211213087,
            -6.06854241,
            -6.246537799)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logn <- ifelse(class != 17, temp[class], log(0.01))
  lowerbound_logn <- min(lowerbound_logn, na.rm = TRUE)
}

#'Estimate manning's n upperbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_upperboundlogn <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(-0.377797559,
            -1.024617813,
            -0.469297566,
            -1.330691863,
            -1.17053266,
            -1.181653683,
            -0.819808414,
            -1.085062606,
            -0.540120016,
            -1.349783005,
            -0.325687448,
            -0.025491812,
            -0.554602291,
            0.501441352,
            1.930473909,
            0.094113848)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logn <- ifelse(class != 17, temp[class], log(0.05))
  upperbound_logn <- max(upperbound_logn, na.rm = TRUE)
}

#'Estimate manning's n SD using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lognSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(1.069101032,
            1.004209047,
            1.051391594,
            1.050420395,
            1.103812115,
            1.050696519,
            1.205687421,
            1.229666517,
            1.101795647,
            1.194241604,
            1.268715794,
            1.260820049,
            1.199816945,
            1.345333569,
            1.724262114,
            1.122652969)

  class <- apply(Wobs, 1, classify_func)
  logn_sd <- ifelse(class != 17, temp[class], 0.761673112)
}

#'Classify river for expert framework
#'
#'@param Wobs observed widths matrix
classify_func <- function(Wobs) {
  lwbar <- mean(log(Wobs), na.rm=TRUE)
  lwsd <- sd(log(Wobs), na.rm= TRUE)

  maxWidth = 6.5
  classes <- c(2.476118144,
               2.864001065,
               3.103015939,
               3.249308032,
               3.284178964,
               3.371669039,
               3.56827873,
               3.664586762,
               3.683922384,
               4.002696788,
               4.031559142,
               4.357733942,
               4.436574004,
               4.921166637,
               5.287893051) #median width of each river type
  index <- ifelse(lwbar > maxWidth, 17, which.min(abs(classes-lwbar))) #17 for big rivers
  index <- ifelse(lwsd >= 0.45, 16, index)  #16 for width-variable rivers, which overrides 'big' rivers
  return(index)
}
