## source simulation file required to gt simulating algorithm
source('rscripts/simulation.R')

#' Triggering Kernel Function
#'
#'  calculates the influence of a given event at the givent time the event is a 2
#' elements list (mi, ti). Inclusive parameter means that if an event is present
#' at time t, then it contributes to the conditional intensity. If inclusive ==
#' F, then events at time t are removed.
#' @param event the event for which we need kernel function
#' @param t time at which we need to calculate the kernel
#' @param K kappa for our social kernel
#' @aparm alpha the exponent for user distribution
#' @param beta the exponent for user influence in social kernel
#' @param mmin the minimum value for powerlaw user distribution
#' @aparma c the cutoff parameter in social kernel
#' @param theta th powerlaw exponent of social kernel
#' @param inclusive whether to consider event at time t for calculation 
kernelFct<- function(event, t, K = 0.024, alpha = 2.016, beta = 0.5, mmin = 1, 
                     c = 0.001, theta = 0.2, inclusive = T) {
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
  return (val)
}


#' Conditional intensity function
#' 
#' The CIF function is necessary to simulate a non-stationary 
#' (NON-HOMOGENUOUS) poisson proces# conditional intensity function
#' in here need to calculate the conditional intensity, by using the kernels
#' @param x the time t till which we need CIF
#' @param history the history of events in cascade
#' @param ... the parameters to be passed to kernelFCT
CIF = function(x, history, ...) {
  subst <- history[history$time <= x,]
  return(sum(kernelFct(event = subst, t = x, ...)))
}

#' Simulate marked Hawkes 
#' 
#' Main function to generate a Hawkes process sequence. It allows intermediary 
#' saves and continuing a stopped simulation. Creates a CSV file with two
#' columns, each row is an event: (magnitude, time)
#' @param K, alpha, beta, mmin, c, theta - parameters of the Hawkes kernel
#' @param M - magnitude of the initial event (in case no initial history 
#'   provided)
#' @param history_init An initial history can be provided (R structure obtained 
#'   from a previous call of this function). Its purpose is to allow custom 
#'   initializations and to continue simulation of stopped processes.
#' @param Tmax - maximum time of simulation.
#' @param filename - file to which save the CSV file with the simulation.
generate_Hawkes_event_series <- function(K = 0.024, alpha = 2.016, beta = 0.5, mmin = 1, c = 0.001, theta = 0.2,
                                         M=10000, Tmax = 10, filename = NULL, history_init = NULL) {
  saveInterval = 1000
  history <- NULL
  
  # initial event, magnitude M and time 0
  if ( !is.null(history_init)) {
    history <- history_init
    t <- history[nrow(history), "time"]
  }
  
  # maybe we have to continue simulation
  if (is.null(filename)) {
    if (is.null(history)) {
      t = 0; 
      history <- as.data.frame(matrix(c(M, t), nrow=1, byrow = TRUE))
      colnames(history) <- c("magnitude", "time")
    }
  } else { 
    if (file.exists(filename)) {
      # means we are continuing the simulation
      # check first if there is some initialized history
      if ( !is.null(history)) {
        warning(sprintf("We gave me both a simulation to continue (yes, file %s exists) and an initial history. Ignoring initial history!", filename))
      }
      history <- read.table(file = filename, header = T)
      colnames(history) <- c("magnitude", "time")
      t <- history[nrow(history), "time"] 
      cat(sprintf("--> Loaded %d events from file, simulation time %.3f.\n", nrow(history), t))
    } else {
      # means we mearly want to save results to file
      # check first if there is some initialized history
      if (is.null(history)) {
        t <- 0
        history <- as.data.frame(matrix(c(M, t), nrow=1, byrow = TRUE))
        colnames(history) <- c("magnitude", "time")
      }
      cat(sprintf("--> Will save progress to history file %s, every %d events!\n", filename, saveInterval))
    }
  }
  
  while(history[nrow(history), "time"] <= Tmax){
    # generate the time of the next event, based on the previous events
    t.lambda.max <- t
    intensityMax <- CIF(x = t.lambda.max, history = history, K, alpha, beta, mmin, c, theta)
    
    # I need this other function, because in the Poisson generation, the time index starts at 0
    # and I need to translate it at the end of my time.
    subsCIF <- function(x, ...) CIF(x + t, ...) 
    # the Tmax = NULL is a dirty trick to convince R to accept my additional parameters
    x <- rnhpoisson(Tmax = NULL, Nmax = 1, LambdaMax = intensityMax, FUN = subsCIF, 
                    history, K, alpha, beta, mmin, c, theta) ## these parameters get transmitted using ...
    t = t + x
    
    # generate the influence of the next event, by sampling the powerlaw distribution of the #retweets
    mag <- generate_user_influence(n = 1, alpha = alpha, mmin = mmin)
    
    # add the next event to the history
    event <- matrix( c(mag, t), nrow=1, byrow = T)
    colnames(event) <- c("magnitude", "time")
    history = rbind(history, event)
    
    cat(sprintf("\rCurrent simulation time: %.3f / %.3f (%d events).", t, Tmax, nrow(history)))
    
    # save progress
    if ((nrow(history) %% saveInterval == 0) && (!is.null(filename))) {
      write.table(x = history, file = filename, sep = "\t", row.names = F, col.names = T)
    }
  }
  
  # if a filename was provided, write the history to the file
  if (!is.null(filename)) {
    write.table(x = history, file = filename, sep = "\t", row.names = F, col.names = T)
  }
  
  cat(sprintf("\n--> Simulation done!\n"))
  return(history)
}

#' Calculate Branching factor
#' 
#' calculates the average number of events generated by an event
#' in our Social kernel
#' useful to detect super-critical regimes
#' @param K kappa for our social kernel
#' @aparm alpha the exponent for user distribution
#' @param beta the exponent for user influence in social kernel
#' @param mmin the minimum value for powerlaw user distribution
#' @aparma c the cutoff parameter in social kernel
#' @param theta th powerlaw exponent of social kernel
#' @return brnaching factor forthe social kernel with parameters passed
getBranchingFactor<- function(K = 0.024, alpha = 2.016, beta = 0.5, mmin = 1, c = 0.001, theta = 0.2) {
  
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

#' Gets total tweets for our Social Kernel with current history
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
    total.tweets = round(dim(history)[1] + a1 / (1 - n.star))
    return (c(total = total.tweets, nstar = n.star, a1 = a1))
  }
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

#' A function to calculate the value of the integral of lambda for
#' Social Kernel
integrateLambda<- function(lower, upper, history, params, mmin = 1) {
  if (length(params) == 4)
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

#'  calculates the negative log-likelihood. 
#'  This is given that the optim function minimizes a function. '
#' Minimizing the -1 * log-likelihood amounts to maximizing the log-likehood. '
neg.log.likelihood <- function(params, history,
                               lowerBound = c(K = 0, beta = 0, c = 0, theta = 0), 
                               upperBound = c(K = Inf, beta = Inf, c = Inf, theta = Inf)
                               ) { 
  
  if (length(params) == 4)
    names(params) <- c("K", "beta", "c", "theta")
  params <- as.list(unlist(params))
  
  ## sanity check
  maxValue <- .Machine$double.xmax - 1
  # if (sum(is.nan(unlist(params)), na.rm = T) > 0) return(maxValue)
  
  ## first check bounds
  # check if we got bounds
  if (is.null(lowerBound)) lowerBound <- rep(-Inf, times = length(params))
  if (is.null(upperBound)) upperBound <- rep(Inf, times = length(params))
  
  ## if all good, continue
  bigT <- max(history$time)
  
  ## get the formula of the log-likelihood * -1.
  ## that in the summation, we eliminated the first event, as there is no background rate in our lambda(t)
  return(integrateLambda(lower = 0, upper = bigT, history = history, params = params) - 
           sum(log(lambda(t = history$time[-1], history = history, params = params))) )
}

#' A function to calculate the dervative in closed form for our Social Kernel
closedGradient <- function(params, history,
                           lowerBound = c(K = 0, beta = 0, c = 0, theta = 0), 
                           upperBound = c(K = Inf, beta = Inf, c = Inf, theta = Inf)){
  if (length(params) == 4)
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

#' A function to specify the constraint.
constraint <- function(params, history,
                       lowerBound = c(K = 0, beta = 0, c = 0, theta = 0), 
                       upperBound = c(K = Inf, beta = Inf, c = Inf, theta = Inf)) {
  if (length(params) == 4)
    names(params) <- c("K", "beta", "c", "theta")
  
  return( log(params[1]) + log(1.1016) - log(1.016 - params[2]) - 
            log(params[4]) - ( params[4] * log(params[3]) ) )
}

#' A function to specify the Jacobian.
jacobian <- function(params, history,
                     lowerBound = c(K = 0, beta = 0, c = 0, theta = 0), 
                     upperBound = c(K = Inf, beta = Inf, c = Inf, theta = Inf)) {
  if (length(params) == 4)
    names(params) <- c("K", "beta", "c", "theta")
  
  return ( c( 1 / params[1], 1 / (1.016 - params[2]), -params[4] / params[3],
              -(1 / params[4]) - log(params[3]) ) )
}

## Functions listed here are helper functions to help us
#' Function to create random and defined starting points
createStartPoints <- function(seed = NULL){
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

#' Function to run ipopt
require(ipoptr) 
fitParameters <- function(start, history){
  
  if (length(start) == 4)
    names(start) <- c("K", "beta", "c", "theta")
  start <- as.list(unlist(start))
  
  opts <-list(print_level = 0, linear_solver = "ma57", max_iter = 10000
              #, "file_print_level"=7, output_file="scripts/ipoptr.out"
  )
 
  constraint_lb <- c(log(.Machine$double.eps))
  constraint_ub <- c(log(1 - .Machine$double.eps))
  
  res <- ipoptr(x0 = c(K = start$K, beta = start$beta, 
                       c = start$c, theta =start$theta),
                eval_f = neg.log.likelihood, 
                eval_grad_f = closedGradient, 
                eval_g= constraint,
                eval_jac_g = jacobian,
                eval_jac_g_structure = list(c(1,2,3,4)),
                lb = c(K = 0, beta = 0, c = 0, theta = 0), 
                ub = c(K = 1, beta = 1.016, c = Inf, theta = Inf),
                constraint_lb = constraint_lb, 
                constraint_ub = constraint_ub,
                opts = opts, # options for algorithm
                history = history,
                lowerBound = c(K = 0, beta = 0, c = 0, theta = 0),
                ## these are parameters for the fn and gr functions
                upperBound = c(K = 1, beta = 1.016, c = Inf, 
                               theta = Inf)
  )
  return (res)
}