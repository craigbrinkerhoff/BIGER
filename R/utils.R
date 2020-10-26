# utility functions

#' Convert coefficient of variation to sigma parameter of lognormal diistribution
#'
#' @param cv Coefficient of variation
#' @export

cv2sigma <- function (cv) {
  sqrt(log(cv^2 + 1))
}

#Functions for classifying rivers----------------------------------------------------------------------

#'Classify river for expert framework
#'
#'@param Wobs observed widths matrix
classify_func <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA
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

#'Classify river for k600 prior assignment
#'
#'@param Wobs observed widths matrix
classify_func_k600 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA
  lwbar <- mean(log(Wobs), na.rm=TRUE)

  classes <- c(log(10),
               log(50),
               log(100))
  index <- ifelse(lwbar < classes[1], 1,
                  ifelse(lwbar >= classes[1] & lwbar < classes[2], 2,
                         ifelse(lwbar >= classes[2] & lwbar < classes[3],3,4)))
  return(index)
}
