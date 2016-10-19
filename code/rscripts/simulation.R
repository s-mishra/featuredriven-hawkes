# Load required libraries
require(poweRlaw)


#' User mark/followers
#' 
#' Samples a power law distribution to genrate user's #followers
#' the exponent alpha=2.016 is derived from a random sample of
#' Twitter data and min value is set to 1.
#' @keywords internal
#' @param n number of marks
#' @param alpha exponent of mark powerlaw
#' @param mmin minimum value to be used in powerlaw
#' @return a list of size n with follwers/mark values
generate_user_influence <- function(n, alpha=2.016, mmin=1){
  user_infl <- conpl$new()
  ## initlaize min and exponent value to generate a powerlaw
  user_infl$setXmin(mmin)
  user_infl$setPars(alpha)
  
  return(dist_rand(user_infl, n))
}

#' Simulate homogeneous poisson process
#' 
#' simulates a homogeneous poisson process with either a maxmium time interval
#' or maximum number of events.
#' algotithm taken from
#' http://transp-or.epfl.ch/courses/OptSim2012/slides/05b-poisson.pdf
#' @param Tmax maximum time for simulation
#' @param Nmax m,aximum number of events to simulate
#' @param lambda the intesity rate for lambda
#' @return a list of all simulated events
rpoisson <- function(Tmax = NULL, Nmax = NULL, lambda) {
  
  # we can have both NULL or both set at the same time.
  if ( !xor(is.null(Tmax), is.null(Nmax)) ) {
    stop("Need to set one (and only one) of Nmax or Tmax")
  }
  
  t = 0
  k = 0
  S = vector()
  
  while (T) {
    r <- runif(1)
    t <- t - log(r) / lambda
    k <- k + 1
    S <- c(S, t)
    
    if (!is.null(Tmax) && (t >= Tmax)) {
      break;
    }
    if (!is.null(Nmax) && (length(S) >= Nmax)) {
      break;
    }
  }
  
  return (S)
}

#' Simulate non-homogeneous poisson process
#' 
#' simulates a non-homogeneous poisson process with either a maxmium time 
#' interval or maximum number of events.
#' algotithm taken from
#' http://transp-or.epfl.ch/courses/OptSim2012/slides/05b-poisson.pdf
#' @param Tmax maximum time for simulation
#' @param Nmax m,aximum number of events to simulate
#' @param FUN funtion to get time dependent intensity at time t
#' @param LamdaMax the intesity rate for rejection sampling
#' @param ... additional parameters required by FUN
#' @return a list of all simulated events
rnhpoisson <- function(Tmax = NULL, Nmax = NULL, FUN, LambdaMax, ...) {
  
  # we can have both NULL or both set at the same time.
  if ( !xor(is.null(Tmax), is.null(Nmax)) ) {
    stop("Need to set one (and only one) of Nmax or Tmax")
  }
  t = 0
  k = 0
  S = vector()
  
  while (T) {
    r <- runif(1)
    t <- t - log(r) / LambdaMax
    
    s <- runif(1)
    thr <- FUN(t, ...)/LambdaMax
    if ( s <= thr) {
      k <- k + 1
      S <- c(S, t)
    }
    
    if (!is.null(Tmax) && (t >= Tmax)) {
      break;
    }
    if (!is.null(Nmax) && (length(S) >= Nmax)) {
      break;
    }
  }
  
  return (S)
}