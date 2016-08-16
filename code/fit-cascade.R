## This script is for running the fitting of Hawkes with IPOPT on a given cascade, you need to give 
## time till the cascsdeto be used for fitting.

# get ipoptr for fitting
require(ipoptr)

#' Function to create random and defined starting points
#' @keywords internal
#' @param seed as the seed if you want to fix the seed for random point generation
create.start.points <- function(seed = NULL){
  set.seed(seed)
  start.K <- runif(n = 5, min = .Machine$double.eps, max = 1)
  start.beta <- runif(n = 5, min = .Machine$double.eps, max = 1.016)
  start.theta <- runif(n = 5, min = .Machine$double.eps, max = 1)
  start.c <- runif(n = 5, min = .Machine$double.eps, max = 1)
  start.K[6:8] <- c(0.1, 0.5, 0.01)
  start.beta[6:8] <- c(0.1, 0.5, 0.01)
  start.theta[6:8] <- c(0.1, 0.5, 0.01)
  start.c[6:8] <- c(0.1, 0.5, 0.01)
  start <- data.frame('K' = start.K, 'beta' = start.beta,
                      'c' = start.c, 'theta' = start.theta)
  return(start)
}

#' A function to calculate the value of the integral of lambda for
#' Power Law Kernel
#' @keywords internal
#' @param lower as the value for time to start integration from
#' @param upper as the value of time as upper limit of integration
#' @param params as parameters for hawkes process
#' @param mmin as the minimum value to be used for magnitude distribution
integrateLambda<- function(lower, upper, history, params, mmin = 1) {
  names(params) <- c("K", "beta", "c", "theta")
  params <- as.list(unlist(params))
  
  ## closed form
  res <- sum(sapply(X = 1:nrow(history), FUN = function(i) {
    val <- (history$magnitude[i] / mmin)^params$beta * 
      (   1 / (params$theta * params$c^params$theta) - 
            1 / (params$theta * (upper + params$c - history$time[i])^params$theta)
      )
    return(val)
  })) * params$K
  
  return(res)
}
kernelFct <- function(event, t, K = 0.024, alpha = 2.016, beta = 0.5, mmin = 1, c = 0.001, theta = 0.2, inclusive = T) {
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]
  
  # f(p_j) part - virality of a video. Constant for a given video
  fun_f <- K
  
  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ beta
  
  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- 1 / (t - ti + c)^(1+theta)
  
  val = fun_f * fun_ro * fun_psi
  val[t<ti] = 0
  val[mi<mmin] = 0
  if (!inclusive) {
    val[t == ti] = 0
    val[mi == mmin] = 0
  }
  
  
  (val)
}
#' Computes the conditional intensity of the non-homogenous Poisson process (
#' lambda(t) ) It is the equivalent of the CIF function in "simulation-model.R",
#' but updated for passing parameters as a list and calculating the intensity
#' for a vector on input (needed for numerical integration).
lambda <-  function(t, history, params = c(K = 0.024, beta = 0.5, c = 0.001, theta = 0.2), inclusive = F, ...) {
  if (length(params) == 4)
    names(params) <- c("K", "beta", "c", "theta")
  params <- as.list(unlist(params))
  
  res <- sapply(X = t, FUN = function(ti) {
    subst <- history[history$time <= ti,]
    return(sum(kernelFct(event = subst, t = ti, 
                         K = params$K, beta = params$beta, c = params$c, theta = params$theta, 
                         inclusive = inclusive, ...)))
  })
  
  return(res)
}

#' Next function calculates the negative
#' log-likelihood. This is given that the optim function minimizes a function. '
#' Minimizing the -1 * log-likelihood amounts to maximizing the log-likehood. '
#' @keywords internal
#' @param params as the parametrs value to use to get likeli-hood
#' @param history as the cascade to be used for evaluation
#' @param lowerBound as the list of values specifying the lower bound on parameters
#' @param upperBound as the list of values specifying the upper bound on parameters
neg.log.likelihood <- function(params, history) { 
  names(params) <- c("K", "beta", "c", "theta")
  params <- as.list(unlist(params))
  
  bigT <- max(history$time)
  
  ## get the formula of the log-likelihood * -1.
  ## that in the summation, we eliminated the first event, as there is no background rate in our lambda(t)
  return(integrateLambda(lower = 0, upper = bigT, history = history, params = params) - 
           sum(log(lambda(t = history$time[-1], history = history, params = params))) )
}

#' constraint function for Power Law Kernel
#' @keywords internal
constraint <- function( params, history) {
  return( log(params[1]) + log(1.1016) - log(1.016 - params[2]) - 
            log(params[4]) - ( params[4] * log(params[3]) ) )
}

#' Jacobian for Power Law Kernel
#' @keywords internal
jacobian <- function (params, history){
  return ( c( 1 / params[1], 1 / (1.016 - params[2]), -params[4] / params[3],
              -(1 / params[4]) - log(params[3]) ) )
}

#' A function to calculate the dervative in closed form for Power Law Kernel.
#' @keywords internal
closed.gradient<- function(params, history){
  names(params) <- c("K", "beta", "c", "theta")
  
  # variables collection
  params <- as.list(unlist(params))
  bigT <- max(history$time)
  n <- dim(history)[1]
  k <- params$K
  beta <- params$beta
  c <- params$c
  theta <- params$theta
  
  # In all the following calculation, res is the part of derivative coming 
  # from bigLambda part and first is derivative coming from 
  # summation of log(lambda(ti))
  
  ### calculating derivative wrt k in closed form
  mi.beta <- history$magnitude ^ beta
  first.part <- 1 / (c ^ theta)
  second.part <- bigT + c - history$time
  second.part.k <- 1 / (second.part ^ theta)
  res <- sum(mi.beta * (first.part - second.part.k)) / theta
  derivK <- ((n-1) / k) - res
  
  ### calculating dervative wrt beta in closed form
  ## MAR: second part is almost identical to the one above, reusing.
  res <- k * sum(mi.beta * log(history$magnitude) * (first.part - second.part.k)) / theta
  
  # we go only from second row as first row has no history before
  # (inner sum goes strictly one less) and we do not have a mu(background rate)
  # numerically (history$magnitude[1:i-1] or (history$time[1:i-1] is zero
  first <- data.frame(t(sapply(X = 2:nrow(history), FUN = function(i) {
    ## MAR: nominator and denominator are basically the same, just nominator has an additional log(m_j)
    denominator <- (((history$time[i]+c) - history$time[1:i-1]) ^-(1+theta)) * (history$magnitude[1:i-1] ^ beta)
    numerator <- denominator * log(history$magnitude[1:i-1])
    return( list(first = sum(numerator) / sum(denominator), denominator = denominator) )
  })))
  
  derivBeta <- sum(unlist(first$first)) - res
  
  ### calculating derivative wrt c in closed form
  first.part <- 1 / (c ^ (1 + theta) )
  second.part.c <- 1 / (second.part ^ (1 + theta) )
  res <- sum(mi.beta * (second.part.c - first.part)) * k
  
  first.c <- sum(sapply(X = 2:nrow(history), FUN = function(i) {
    
    denominator <- first$denominator[[i-1]]
    numerator <- sum((((history$time[i]+c) - history$time[1:i-1]) ^-(2+theta))
                     * -(1+theta) * (history$magnitude[1:i-1] ^ beta))
    
    # denominator <- sum(((history$time[i]+c) - history$time[1:i-1]) ^-(1+theta) 
    #                    * (history$magnitude[1:i-1] ^ beta))
    
    return(numerator / sum(denominator) )
  }))
  derivC <- first.c - res
  
  #### calculating derivative wrt theta in closed form
  first.part <- second.part ^ (-theta) * (theta * log(second.part) + 1)
  second.part.theta <- c ^ (-theta) * ((theta * log(c)) + 1)
  res1 <- sum(mi.beta * (first.part - second.part.theta)) * k / (theta ^ 2)
  
  res <- sum(sapply(X = 1:nrow(history), FUN = function(i) {
    
    # calculating in three parts so that can be debugged easily and 
    # expressed in terms of smaller terms
    
    mi.beta <- (history$magnitude[i]) ^ beta 
    
    first.part <- ((bigT + c) - history$time[i]) ^ (-theta) * 
      ((theta * log((bigT + c) - history$time[i])) + 1)
    
    second.part <- c ^ (-theta) * ((theta * log(c)) + 1)
    
    return(mi.beta * (first.part - second.part))
  })) * k / (theta ^ 2)
  
  first.theta <- -sum(sapply(X = 2:nrow(history), FUN = function(i) {
    
    denominator <- first$denominator[[i-1]]    
    numerator <- log(history$time[i] + c - history$time[1:i-1]) * denominator
    return( sum(numerator) / sum(denominator) )
  }))
  derivTheta <- first.theta - res
  
  ## added a -1 because we are minimizing log-likelihood
  return(-1 * c(derivK, derivBeta, derivC, derivTheta))
}

#' Function that fit's a hawkes process to a given cascade
#' @param history as the cascade to be fitted
#' @param starting intialization point 
fit.hawkes <- function(history, start){
  ## run the fitting
  tryCatch({
    # call ipopt
    # get starting points for intializations
    names(start) <- c("K", "beta", "c", "theta")
    start <- as.list(unlist(start))
    
    constraint_lb <- c(log(.Machine$double.eps))
    constraint_ub <- c(log(1 - .Machine$double.eps))
    
    opts <-list(print_level = 0, linear_solver = "ma57", max_iter = 10000)
    
    res <- ipoptr(x0 = c(K = start$K, beta = start$beta, 
                         c = start$c, theta =start$theta),
                  eval_f = neg.log.likelihood, 
                  eval_grad_f = closed.gradient, 
                  eval_g= constraint,
                  eval_jac_g = jacobian,
                  eval_jac_g_structure = list(c(1,2,3,4)),
                  lb = c(K = 0, beta = 0, c = 0, theta = 0), 
                  ub = c(K = 1, beta = 1.016, c = Inf, theta = Inf),
                  constraint_lb = constraint_lb, 
                  constraint_ub = constraint_ub,
                  opts = opts, # options for algorithm
                  history = history
    )
    ## when here, we have finished running IPOPT
    retval <- list(solver.output = res, obj.val = res$objective, params = res$solution)
    return(retval)
  }, error = function(err) {
    retval <- list(obj.val = NA, params = NA, solver.output = NA)
    return(retval)
  })
  
}
#' the get_n function calculates the average number of events generated by an event
#' in Power Law Kernel
#' useful to detect super-critical regimes
getBranchingFactor <- function(K = 0.024, alpha = 2.016, beta = 0.5, mmin = 1, c = 0.001, theta = 0.2) {
  
  if (beta >= (alpha - 1) ) {
    warning("The closed expression calculated by this function does NOT hold for beta >= alpha - 1")
    
    return (Inf)
  }
  
  if (theta <= 0 ) {
    warning(sprintf("The closed expression calculated by this function does NOT hold for theta <= 0 (K=%.4f, beta=%.2f, theta=%.2f)", K, beta, theta))
    
    return (Inf)
  }
  
  n0 = (K * (alpha - 1)) / (alpha - 1 - beta)
  int = 1 / ( theta * c^theta)
  n = n0 * int
  
  return (n) 
}

#' Gets total tweets for a power law kernel with current history
#' 
getTotalEvents<- function(history, bigT, K = 0.024, alpha = 2.016, beta = 0.5, mmin = 1, c = 0.001, 
                               theta = 0.2 ){
  n.star <- getBranchingFactor(K, alpha, beta, mmin, c, theta)
  if (n.star >= 1) {
    warning(sprintf("Branching Factor greater than 1, not possible to predict the size(super critical)"))
    
    return (c(total = Inf, nstar = n.star, a1 = 'NA'))
    
  } else{
    # getting T from data -- WRONG!!! bigT is the observed time 
    # bigT <- max(history$time)
    # calculating the expected size of first level of descendants
    a1 <-  sum(sapply(X = 1:nrow(history), FUN = function(i) {
      val <- (history$magnitude[i]) ^ beta / (theta * ((bigT + c - history$time[i]) ^ theta))
      return(val)
    })) * K
    # calculating the final value
    total.tweets = dim(history)[1] + a1 / (1 - n.star)
    return (c(total = total.tweets, nstar = n.star, a1 = a1))
  }
}

#' function used to get prediction for a set of fitted values
#' @param tweets as the cascade to be predicted
#' @param pred.time as the time after which tomake prediction
get.prediction <- function(tweets, pred.time){
  colnames(tweets) <- c('magnitude', 'time')
  history <- tweets[with(tweets, time <= pred.time),]
  start <- create.start.points()
  mle <- list()
  for (s in 1:nrow(start)) {
    myrow <- list(obj.val = NA, params = NA )
    retval <- fit.hawkes(history, start[s,])
    myrow <- retval
      if (retval$solver.output$status != 0 & 
            getBranchingFactor(K = retval$params[1], beta = retval$params[2], c = retval$params[3], theta = retval$params[4]) >= 1) {
        myrow$obj.val <- NA
      }
  
    mle <- c(mle, list(myrow))
  }
  
  mle <- data.frame(t(sapply(mle, FUN = function(el) return(el))))
  
  # try as if we produced all infinite or NA as mle which is because 
  # of not possible values then we will get error here and escape
  tryCatch({
    # which.min works on NA by default!!
    optim.res <- unlist(mle$params[which.min(mle$obj.val)])
    pred <- getTotalEvents(history, bigT = pred.time,
                           K = optim.res[1], alpha = 2.016, 
                           beta = optim.res[2], mmin = 1, 
                           c = optim.res[3], theta = optim.res[4])
    return(pred)
    
  }, error = function(err) {
    return(NA)
  })
  
}

### to check providing a test data
###load(file = 'test.Rdata')
### here 300 is time for which we want to train
### pred <- get.prediction(test,300)


