######################
#Functions for default prior setup for algorithm (n and A0 functions from Brinkerhoff etal 2020 --> geoBAM-Expert algorithm)
######################


#k prior calculations--------------------------------------------------------------------------------------------------
#' logk_hat prior estimation
#' 
#' Estimate k_hat prior hyperparameter using swot data
#'
#' @param Sobs Observed S,as a space-down, time-across matrix
#' @param priorQ mean annual flow estimate, single number
#' @export
estimate_logk <- function(Sobs, priorQ){
  Sobs[Sobs <=0] <- NA

  slope <- mean(Sobs, na.rm=T)
  
  #prior via Brinkerhoff etal 2019 rating curves and Qwbm prior
  khat <- rep(log(62.82*(9.8*slope)^(7/16)*(0.252*priorQ^0.388)^(9/16)*(0.276*priorQ^0.164)^(1/4)), ncol(Sobs))
}

#' logk_sd prior estimation
#' 
#' Estimate k sd prior hyperparameter using swot data
#'
#' @param Sobs Observed S,as a space-down, time-across matrix
#' @export
estimate_logksd <- function(Sobs){
  Sobs[Sobs <=0] <- NA

  #ksd <- rep(log(3.45), ncol(Sobs)) #standard error for CHAINSAW/ES MODEL
  ksd <- rep(0.30, ncol(Sobs)) #standard error for CHAINSAW/ED MODEL
}

# Prior calculation using geoBAM-Expert classification framework------------------------------------------------------------------
#class 17 are 'big' rivers
#class 16 are highly width variable
#These are now using a lookup table generated from additional filtering for the geobam prior dataset FYI. See ~/ngoing_Projects/geoBAM_update_Summer_2021 for the script and lookup tables

#' logA0_hat prior estimation
#' 
#' Estimate median cross-sectional area prior hyperparameter using swot data
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
            4.736198448,
            4.990432587,
            4.844063071,
            5.369441357,
            5.422742508,
            5.523458921,
            5.774532256,
            6.446836611,
            6.527953643,
            6.873163834,
            7.13630945,
            7.996317232,
            8.743785303,
            4.428426599)

  class <- apply(Wobs, 1, classify_func)
  logA0_hat <- ifelse(class != 17, temp[class], -0.2918 + 1.6930 * lwbar - 0.1887 * lwsd) #repeat for each sptial unit
  #global r2: 0.907
}

#' lowerbound_A0 prior estimation
#' 
#' Estimate median cross-sectional area lowerbound prior hyperparameter using swot data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification (p5)
  temp <- c(1.627714731,
            2.437710905,
            2.643581022,
            2.944438979,
            2.862725217,
            3.692086845,
            3.890215272,
            3.46067656,
            3.613981618,
            4.042002548,
            4.063218741,
            4.203662669,
            5.55624353,
            5.688828434,
            5.928895938,
            1.942398159)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_A0 <- ifelse(class != 17, exp(temp[class]), exp(5.928895938)) #class 15 p5 value
  lowerbound_A0 <- min(lowerbound_A0, na.rm = TRUE)
}

#' Upperbound_A0 prior estimation
#' 
#' Estimate median cross-sectional area upperbound prior hyperparameter using swot data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification (p95)
  temp <- c(6.592711677,
            7.007838922,
            6.800543165,
            6.782192056,
            6.912893281,
            7.882304158,
            7.461640392,
            7.800850866,
            7.288859014,
            8.178015738,
            8.139949251,
            8.289026216,
            8.585211456,
            9.178526329,
            10.48217535,
            8.994142514)

  class <- apply(Wobs, 1, classify_func)
  upperbound_A0 <- ifelse(class != 17, exp(temp[class]), exp(10.48217535)) #class 15 p95 value
  upperbound <- max(upperbound_A0, na.rm = TRUE)
}

#' logA0_sd prior estimation
#' 
#' Estimate median cross-sectional area SD prior hyperparameter using swot data
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
            1.365242448,
            1.258756782,
            1.302566576,
            1.311046184,
            1.190542029,
            1.337220271,
            1.21047899,
            1.359096608,
            1.245863462,
            1.287297345,
            1.109424148,
            1.144803558,
            1.454662355,
            2.233590129)

  class <- apply(Wobs, 1, classify_func)
  logA0_sd <- ifelse(class != 17, temp[class], 0.58987527)
}

#' logn_hat prior estimation
#' 
#' Estimate manning's n prior hyperparameter using swot data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix
#' @param Sobs Observed s, as a space-down, time-across matrix
#' @export
estimate_logn <- function(Wobs, Sobs) {
  Sobs[Sobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lsbar <- apply(log(Sobs), 1, mean, na.rm = TRUE)
  #Expert classification
  temp <- c(-2.985277801,
            -3.218074003,
            -3.085265829,
            -3.316850463,
            -3.318215717,
            -3.143583451,
            -3.456243594,
            -3.545592423,
            -3.240085716,
            -3.401877538,
            -3.278677422,
            -3.372948626,
            -3.395070774,
            -3.30626809,
            -3.455433411,
            -3.318541261)

  class <- apply(Wobs, 1, classify_func)
  logn_hat <- ifelse(class != 17, temp[class], -0.1636 + 0.4077 * lsbar) #repeat for each spatial unit
  #Global r2: 0.631
}

#' logn_sd prior estimation
#' 
#' Estimate manning's n SD prior hyperparameter using swot data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lognSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(1.029408933,
            1.004209047,
            1.047755392,
            1.055423371,
            1.091631839,
            1.059532393,
            1.205687421,
            1.229666517,
            1.101795647,
            1.194241604,
            1.225854067,
            1.167193097,
            1.152000583,
            1.285632655,
            1.343539617,
            1.102026911)

  class <- apply(Wobs, 1, classify_func)
  logn_sd <- ifelse(class != 17, temp[class], 0.761673112) #standard error from regression in estimate_logn
}
