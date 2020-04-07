#' Group WQS Regression
#'
#' This function fits a group weighted quantile sum (GWQS) regression model.
#'
#' @param y A vector containing outcomes for validation.
#' @param y.train A vector containing outcomes for training. If left as NULL the validation data will be used for training as well.
#' @param x A matrix of component data for validation.
#' @param x.train A matrix of component data for training. If left as NULL the validation data will be used for training as well.
#' @param z A vector or matrix of covariates for validation.
#' @param z.train A vector or matrix of covariates for training. If left as NULL the validation data will be used for training as well.
#' @param x.s A vector of the number of components in each index.
#' @param B The number of bootstrap samples, must be 1 or more.
#' @param n.quantiles The number of quantiles to apply to data.
#' @param pars A vector of initial values.
#' @param func The objective function to be used (must match outcome data type); currently only fun args "continuous" or "binary" are supported.
#' @param ineqLB Vector of lower bounds for betas and weights.
#' @param ineqUB Vector of upper bounds for betas and weights.
#' @param tol Tolerance level for bootstrap convergence.
#' @param delta Step size for bootstrap procedure.
#'
#' @return A list of 2 containing the GWQS estimate based on calculated weights
#' and the GWQS model fit to validation data
#'
#' @examples
#' data("WQSdata")
#' group_list <- list(c("X1", "X2", "X3"), c("X4", "X7"), c("X5", "X6", "X9", "X8"))
#' x.s <- make.x.s(WQSdata, 3, group_list)
#' X <- make.X(WQSdata, 3, group_list)
#' Y <- WQSdata$y
#' results <- wqs.fit(y = Y, x = X, x.s = x.s, B=1, func = "continuous")
#'
#' @import stats
#' @import Rsolnp
#' @import glm2
#'
#' @export
wqs.fit <- function(y, y.train = NULL, x, x.train = NULL, z = NULL, z.train = NULL, x.s, B=100, n.quantiles=4, pars = NULL, func,
                    ineqLB = NULL,
                    ineqUB = NULL,
                    tol = 1e-6,
                    delta = 1e-6){

  if (!is.null(y.train)){
      # Check training data
      check_train <- check_xyz(x.train, y.train, z.train)
      if(check_train[1])stop("x.train must be a matrix")
      if(check_train[2])stop("check dimensions of training data")
      if(check_train[3])stop("check dimensions of training data")
 }

  # Check validation data
  # check_valid <- check_xyz(x, y, z)
  # if(check_valid[1])stop("x must be a matrix")
  # if(check_valid[2])stop("check dimensions of validation data")
  # if(check_valid[3])stop("check dimensions of validation data")

  # Other checks
  if(n.quantiles < 2 | n.quantiles > 10)stop("n.quantiles must be at least 2 and no greater than 10")
  func_check <- ifelse(func == "continuous" | func == "binary", 0, 1)
  if(func_check)stop("Only 'continuous' and 'binary' are currently supported")


  K <- length(x.s) # Total number of component groups
  C <- dim(x)[2] # Total number of components

  ### Defining pars, function, and all training variables if left as NULL, for Z and no Z cases

  if (is.null(z)){         # No Z case

    if (is.null(pars)){    # Defining the vector of initial values if pars left as NULL
      w_inits <- numeric()
      for (i in 1:K){
        w_inits <- c(w_inits, rep(1/x.s[i], x.s[i]))
      }
      pars <- c(rep(0, K+1), w_inits)
    }

    if (length(pars) != sum(K+1, C)){ stop("The length of pars much equal the total number of betas and weights")}


    if (is.null(y.train)){   # Setting all data as training if _.train variables left as NULL
      y.train <- y
      x.train <- x
    }

    if (func == "continuous"){     # Setting Objective Function
      fun <- cont_z
    } else if (func == "binary"){
      fun <- bin_z
    }

  } else {                # Z case

    p <- dim(z)[2] # Total number of covariates

    if (is.null(pars)){        # Defining the vector of initial values if pars left as NULL
      w_inits <- numeric()
      for (i in 1:K){
        w_inits <- c(w_inits, rep(1/x.s[i], x.s[i]))
      }
      pars <- c(rep(0, K+1), w_inits, rep(.1, p))
    }

    if (length(pars) != sum(K+1, C, p)){ stop("The length of pars much equal the total number of betas, weights, and covariates")}

    if (is.null(y.train)){   # Setting all data as training if _.train variables left as NULL
      y.train <- y
      x.train <- x
      z.train <- z
    }

    if (func == "continuous"){    # Setting Objective Function

      fun <- cont
    } else if (func == "binary"){

      fun <- bin
    }

  }
  ### Defining eqB variable
    eqB <- rep(1, K)
  ### Defining ineqLB variable if left as NULL
  if (is.null(ineqLB)){
    ineqLB <- rep(-2, K)
  }
  if(length(ineqLB) != K){ stop("The length of ineqLB must equal the number of component groups")}
  ineqLB <- c(ineqLB, rep(0,C))
  ### Defining ineqUB variable if left as NULL
  if (is.null(ineqUB)){
    ineqUB <- rep(2, K)
  }
  if(length(ineqUB) != K){ stop("The length of ineqUB must equal the number of component groups")}
  ineqUB <- c(ineqUB, rep(1,C))
  ### Taking quantiles
  q <- quantile.fn(x, n.quantiles)
  q.train <- quantile.fn(x.train, n.quantiles)

  ###
  if (B > 1) {

    if (is.null(z)){

      result <- matrix(0, nrow = B, ncol = length(pars))
      colnames(result) <- names(pars)
      convergence <- rep(0, B) #0 indicates convergence; 1 or 2 indicates non-convergence

      betas <- numeric()
      tests <- numeric()
      pvals <- numeric()

      #---------------------------- BOOTSTRAP ROUTINE -----------------------------------
      for (b in 1:B) {
        # draw random sample (of same size as training data set) with replacement
        samp <- sample(1:length(y.train), replace = TRUE)

        y.b <- as.vector(y.train[samp])
        q.b <- q.train[samp,]

        result.b <- solnp(pars, fun, eqfun = lincon1_z, eqB, ineqfun = ineq_z, ineqLB, ineqUB, LB=NULL, UB=NULL, q.b, y.b, x.s,
                          control = list(tol = tol,delta = delta, trace = 0))

        result[b,] <- result.b$pars
        convergence[b] <- result.b$convergence

        w <- result.b$pars[(K+2):(K+1 + sum(x.s))]
        fit <- wqs.fit.internal_z(q.b, y.b, w, x.s, func)$fit
        b.b <- fit$coefficients[-1]
        t.b <- summary(fit)$coefficients[-1,3]
        p.b <- summary(fit)$coefficients[-1,4]
        betas <- rbind(betas, b.b)
        tests <- rbind(tests, t.b)
        pvals <- rbind(pvals, p.b)
      }

      out1 <- list(result, convergence, betas, tests, pvals)

      names(out1) <- c("params", "convergence", "betas", "tests", "pvals")


      wts.matrix <- out1$params[, (K+2):(K+1+ sum(x.s))]
      test_stat <- out1$tests

      weights <- teststat.fn2(wts.matrix, test_stat, x.s)

      final <- wqs.fit.internal_z(q, y, weights, x.s, func)

      return(final)

    } else {


      # initialize matrix for parameter estimates (from estimation step)
      result <- matrix(0, nrow = B, ncol = length(pars))
      colnames(result) <- names(pars)
      convergence <- rep(0, B) #0 indicates convergence; 1 or 2 indicates non-convergence

      betas <- numeric()
      tests <- numeric()
      pvals <- numeric()

      #---------------------------- BOOTSTRAP ROUTINE -----------------------------------
      for (b in 1:B) {
        # draw random sample (of same size as training data set) with replacement
        samp <- sample(1:length(y.train), replace = TRUE)

        y.b <- as.vector(y.train[samp])
        q.b <- q.train[samp,]
        z.b <- as.matrix(z.train[samp,])
        rownames(z.b) <- NULL

        result.b <- solnp(pars, fun, eqfun = lincon1, eqB, ineqfun = ineq, ineqLB, ineqUB, q.b, z.b, y.b, x.s, LB = NULL, UB = NULL,
                          control = list(tol = tol, delta = delta, trace = 0))

        result[b,] <- result.b$pars
        convergence[b] <- result.b$convergence

        w <- result.b$pars[(K+2):(K+1 + sum(x.s))] # w: estimated weights from training data (boot strap samples)
        fit <- wqs.fit.internal(q.b, z.b, y.b, w, x.s, func)$fit
        # Remove first and last coefficients for intercept and Z - assume only 1 Z for now
        b.b <- fit$coefficients[-c(1,(K+2))]
        t.b <- summary(fit)$coefficients[-c(1,(K+2)),3]
        p.b <- summary(fit)$coefficients[-c(1,(K+2)),4]
        betas <- rbind(betas, b.b)
        tests <- rbind(tests, t.b)
        pvals <- rbind(pvals, p.b)
      }

      out1 <- list(result, convergence, betas, tests, pvals)


      names(out1) <- c("params", "convergence", "betas", "tests", "pvals")
      #return(out)

      wts.matrix <- out1$params[, (K+2):(K+1+ sum(x.s))]
      test_stat <- out1$tests

      weights <- teststat.fn2(wts.matrix, test_stat, x.s)


      final <- wqs.fit.internal(q, z, y, weights, x.s, func)

      return(final)
    }

  } else if (B == 1) {

    if (is.null(z)){

      result <- solnp(pars, fun, eqfun = lincon1_z, eqB,
                      ineqfun = ineq_z, ineqLB, ineqUB, LB = NULL, UB = NULL, q = q.train, y = y.train, x.s = x.s,
                      control = list(tol = tol, delta = delta, trace = 0))

      par.estimates <- result$pars

      w_names <- character()
      for (i in 1:K){
        w_names <- c(w_names, paste0("w", i, 1:x.s[i]))
      }
      #names(par.estimates) <- c("b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", paste0("w1", 1:x.s[1]), paste0("w2", 1:x.s[2]), paste0("w3", 1:x.s[3]), paste0("w4", 1:x.s[4]), paste0("w5", 1:x.s[5]), paste0("w6", 1:x.s[6]), paste0("w7", 1:x.s[7]), paste0("z", 1:p))
      names(par.estimates) <- c(paste0("b", 0:K), w_names)
      round(par.estimates, 2)
      weights <- par.estimates[(K+2):(K+1+sum(x.s))]

      final <- wqs.fit.internal_z(q, y, weights, x.s, func)
      return(final)

    } else {

      result <- solnp(pars, fun, eqfun = lincon1, eqB,
                      ineqfun = ineq, ineqLB, ineqUB, LB = NULL, UB = NULL, q = q.train, y = y.train, x.s = x.s,
                      z = z.train, control = list(tol = tol, delta = delta, trace = 0))

      par.estimates <- result$pars

      w_names <- character()
      for (i in 1:K){
        w_names <- c(w_names, paste0("w", i, 1:x.s[i]))
      }
      #names(par.estimates) <- c("b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", paste0("w1", 1:x.s[1]), paste0("w2", 1:x.s[2]), paste0("w3", 1:x.s[3]), paste0("w4", 1:x.s[4]), paste0("w5", 1:x.s[5]), paste0("w6", 1:x.s[6]), paste0("w7", 1:x.s[7]), paste0("z", 1:p))
      names(par.estimates) <- c(paste0("b", 0:K), w_names, paste0("z", 1:p))
      round(par.estimates, 2)
      weights <- par.estimates[(K+2):(K+1+sum(x.s))]

      final <- wqs.fit.internal(q, z, y, weights, x.s, func)
      return(final)
    }
  } else {

    stop("B parameter must have positive, nonzero value")

  }
}

###########################################################################

########## Helper Functions ##########
#----------------------------------------------------------------------------------------------
# Function to perform preliminary check of the data
#----------------------------------------------------------------------------------------------
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
check_xyz <- function(x,y,z){
  class_x <- ifelse(class(x)== "matrix" | class(x)== "data.frame", 0, 1)

  n <- length(y)
  n.x <- dim(x)[1]
  dim_x <- ifelse(n != n.x, 1, 0)

  if(is.null(z)){dim_z <- 0}
  else{
    if(class(z) == "matrix") n.z <- dim(z)[1]
    if(class(z) == "vector") n.z <- length(z)[1]
    dim_z <- ifelse(n != n.z, 1, 0)
  }

  return(c(class_x, dim_x, dim_z))
}

##### Continuous Objective Function
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
cont <- function(param, q, z, y, x.s){
  Large <- 1000000000
  K <- length(x.s)
  b0 <- param[1] 	# intercept


  i <- 1
  w.start <- K+2
  w.end <- K+1+x.s[i]
  q.start <- i        # Starting point is different for q
  q.end <- x.s[i]
  mod <- "b0"
  while(i <= K){
    assign(paste0("b", i), param[(i+1)])
    assign(paste0("w", i), param[w.start:w.end])
    assign(paste0("q", i), q[,q.start:q.end])
    mod <- paste(mod, "+", paste0("b", i),"*", paste0("q", i), "%*%", paste0("w", i), sep="")
    i <- i + 1
    w.start <- w.end + 1
    w.end <- w.end + x.s[i]
    q.start <- q.end + 1
    q.end <- q.end + x.s[i]
  }

  #p <- dim(z)[2]	# of covariates
  #theta <- param[(4+C1+C2):(3+C1+C2+p)] 	# parameters for covariates (length p)

  # Add covariate to model
  C <- length(param)
  theta <- param[C]   # assume only one covariate for now
  mod <- paste(mod, "+", "z", "%*%", "theta", sep="")

  ls <- numeric() # initialize space

  #mu <- b0 + b1*q1%*%w1 + b2*q2%*%w2+ z%*%theta

  mu <- eval(parse(text=mod))
  ls <- (y-mu)**2
  leastsq <- sum(ls)
  # Check for NaN
  return(ifelse(is.nan(leastsq), Large, leastsq))
  #return(leastsq) # minimizes objective fn and we want to minimize least squares
}

##### Binary Objective Function
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
bin <- function(param, q, z, y, x.s){
  Large <- 1000000000
  K <- length(x.s)
  b0 <- param[1] 	# intercept

  i <- 1
  w.start <- K+2
  w.end <- K+1+x.s[i]
  q.start <- i        # Starting point is different for q
  q.end <- x.s[i]
  mod <- "b0"
  while(i <= K){
    assign(paste0("b", i), param[(i+1)])
    assign(paste0("w", i), param[w.start:w.end])
    assign(paste0("q", i), q[,q.start:q.end])
    mod <- paste(mod, "+", paste0("b", i),"*", paste0("q", i), "%*%", paste0("w", i), sep="")
    i <- i + 1
    w.start <- w.end + 1
    w.end <- w.end + x.s[i]
    q.start <- q.end + 1
    q.end <- q.end + x.s[i]
  }

  # Add covariate to model

  theta <- param[(sum(x.s)+K+2):(length(param))]
  term <- paste(mod, "+", "z", "%*%", "theta", sep="")

  logl <- numeric() # initialize space
  eterm <- eval(parse(text=term))
  mu <- 1/(1+exp(-eterm))
  logl <- t(y)%*%log(mu) + t(1 - y)%*%log(1-mu)
  loglik <- sum(logl)
  # Check for NaN
  return(ifelse(is.nan(loglik), Large, -loglik))
  #return(-loglik)  # package minimizes objective function so we return neg. log-likelihood
}

##### Linear constraint (allows us to constrain weights to sum to 1)
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
lincon1 <- function(param, q, z, y, x.s){
  K <- length(x.s)
  i <- 1
  start <- K+2
  end <- K+1+x.s[i]
  ret <- rep(0, K)
  while(i <= K){
    weights <- param[start:end]
    assign(paste0("sum", i), sum(weights))
    ret[i] <- eval(parse(text=paste0("sum", i)))
    i <- i + 1
    start <- end + 1
    end <- end + x.s[i]
  }
  #return(c(sum1, sum2, sum3))
  return(ret)
}

##### Inequality constraints (allows us to put constraints on b1 and the weights)
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
ineq <- function(param, q, z, y, x.s){
  K <- length(x.s)
  i <- 1
  start <- K+2
  end <- K+1+x.s[i]
  ret.b <- rep(0, K)
  ret.w <- numeric()
  while(i <= K){
    assign(paste0("b", i), param[(i+1)])
    assign(paste0("weights", i), param[start:end])
    ret.b[i] <- eval(parse(text=paste0("b", i)))
    ret.w <- c(ret.w, eval(parse(text=paste0("weights", i))))
    i <- i + 1
    start <- end + 1
    end <- end + x.s[i]
  }
  #return(c(b1, b2, b3, weights1, weights2, weights3))
  return(c(ret.b, ret.w))
}

###### Test Stat -  Calculate weights based on relative test statistic
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
teststat.fn2 <- function(wts, test_stat, x.s){
  Sb = abs(test_stat)  # take absolute value so we can calculate relative strength of test statistic
  B <- dim(Sb)[1]
  Sb.s <- apply(Sb, 2, sum)
  Sb.m <- matrix(rep(Sb.s,each=B),nrow=B)
  signal <- Sb/Sb.m
  w <- vector()

  # Apply weights by index
  K <- length(x.s)
  i <- 1
  start <- 1
  end <- x.s[i]
  while(i <= K){
    #sig <- mean(signal[,i])
    #w[start:end]  <- colSums(signal[,i]*wts[,start:end])
    w[start:end]  <- signal[,i] %*% wts[,start:end]
    i <- i + 1
    start <- end + 1
    end <- end + x.s[i]
  }
  #for(i in 1:K){
  #   wts[,start:end]
  #}
  #sig <- apply(signal, 1, mean)
  #w <- colSums(sig*wts)
  return(w)
}

##### WQS Fit Internal (Fits weights from individual bootstraps and does final fit)
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
wqs.fit.internal <- function(q, z, y, w, x.s, func){
  K <- length(x.s)
  i <- 1
  start <- 1
  end <- x.s[i]
  #mod <- "b0"
  WQS <- numeric()
  temp.n <- "y"
  while(i <= K){
    assign(paste0("w", i), w[start:end])
    assign(paste0("q", i), q[,start:end])
    assign(paste0("WQS", i), eval(parse(text=paste0("q", i))) %*% eval(parse(text=paste0("w", i))))
    WQS <- cbind(WQS, eval(parse(text=paste0("WQS", i))))
    temp.n <- c(temp.n, paste0("WQS", i))
    #mod <- paste(mod, "+", paste0("b", i),"*", paste0("q", i), "%*%", paste0("w", i), sep="")
    i <- i + 1
    start <- end + 1
    end <- end + x.s[i]
  }
  temp <- data.frame(cbind(y, WQS, z))
  temp.n <- c(temp.n, paste0("Z", seq(1, dim(z)[2],1)))
  names(temp) <- temp.n
  if (func == "continuous"){
    fit <- glm2(y ~ ., data = temp, family = "gaussian"(link = identity))
  } else if (func == "binary") {
    fit <- glm2(y ~ ., data = temp, family = "binomial"(link = logit))
  }
  #fit <- lm(y ~ ., data = temp)

  out <- list(temp[,-1], fit)
  names(out) <- c("WQS", "fit")
  return(out)
}

##### Quantile function
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
quantile.fn <- function(data, n.quantiles){
  q <- matrix(0, dim(data)[1], dim(data)[2])
  I <- dim(data)[2]
  for(i in 1:I){
    a <- rank(data[,i], ties.method = "first")
    q[,i] <- cut(a, quantile(a, probs = c(0:n.quantiles/n.quantiles)), include.lowest=TRUE)}
  q <- q-1
  colnames(q) <- colnames(data)
  return(q)
}

######## Helper Functions (No Covariate versions) ########

##### Continuous Objective Function (no Z)
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
cont_z <- function(param, q, y, x.s){
  Large <- 1000000000
  K <- length(x.s)
  b0 <- param[1] 	# intercept
  #C <- length(param)
  #w <- param[(K+2):C]
  #b <- param[1:(K+1)]

  i <- 1
  w.start <- K+2
  w.end <- K+1+x.s[i]
  q.start <- i        # Starting point is different for q
  q.end <- x.s[i]
  mod <- "b0"
  while(i <= K){
    assign(paste0("b", i), param[(i+1)])
    assign(paste0("w", i), param[w.start:w.end])
    assign(paste0("q", i), q[,q.start:q.end])
    mod <- paste(mod, "+", paste0("b", i),"*", paste0("q", i), "%*%", paste0("w", i), sep="")
    i <- i + 1
    w.start <- w.end + 1
    w.end <- w.end + x.s[i]
    q.start <- q.end + 1
    q.end <- q.end + x.s[i]
  }

  #p <- dim(z)[2]	# of covariates
  #theta <- param[(4+C1+C2):(3+C1+C2+p)] 	# parameters for covariates (length p)

  ls <- numeric() # initialize space

  #mu <- b0 + b1*q1%*%w1 + b2*q2%*%w2+ z%*%theta
  #mu <- b0 + b1*q1%*%w1 + b2*q2%*%w2 + b3*q3%*%w3
  #mu <- b0 + b*q%*%w
  mu <- eval(parse(text=mod))
  ls <- (y-mu)**2
  leastsq <- sum(ls)
  return(ifelse(is.nan(leastsq), Large, leastsq)) # Check for NaN, returns least square estimate
}

##### Binary Objective Function (No Z)

# Old Version
# bin_z <- function(param, q, y, x.s){
#   Large <- 1000000000
#   c <- dim(q)[2]	# number of components
#   #p <- dim(z)[2]	# number of covariates
#   b0 <- param[1] 	# intercept term
#   b1 <- param[2] 	# coefficient for WQS term
#   w <- param[3:(2+c)] 			# vector of weights (length c)
#   #theta <- param[(3+c):(2+c+p)] 	# parameters for covariates (length p)
#   logl <- numeric() # initialize space
#   term <- b0 + b1*q%*%w# + z%*%theta
#   mu <- 1/(1+exp(-term))
#   logl <- t(y)%*%log(mu) + t(1 - y)%*%log(1-mu)
#   loglik <- sum(logl)
#   # Check for NaN
#   return(ifelse(is.nan(loglik), Large, -loglik)) # package minimizes objective function so we return neg. log-likelihood
# }

#' title
#'
#' description
#'
#' deets
#'
#' @noRd
bin_z <- function(param, q, y, x.s){
  Large <- 1000000000
  K <- length(x.s)
  b0 <- param[1] 	# intercept

  i <- 1
  w.start <- K+2
  w.end <- K+1+x.s[i]
  q.start <- i        # Starting point is different for q
  q.end <- x.s[i]
  mod <- "b0"
  while(i <= K){
    assign(paste0("b", i), param[(i+1)])
    assign(paste0("w", i), param[w.start:w.end])
    assign(paste0("q", i), q[,q.start:q.end])
    mod <- paste(mod, "+", paste0("b", i),"*", paste0("q", i), "%*%", paste0("w", i), sep="")
    i <- i + 1
    w.start <- w.end + 1
    w.end <- w.end + x.s[i]
    q.start <- q.end + 1
    q.end <- q.end + x.s[i]
  }

  # Add covariate to model
  #theta <- param[(sum(x.s)+K+2):(length(param))]
  #term <- paste(mod, "+", "z", "%*%", "theta", sep="")

  logl <- numeric() # initialize space
  eterm <- eval(parse(text=mod))
  #eterm <- eval(parse(text=term))
  mu <- 1/(1+exp(-eterm))
  logl <- t(y)%*%log(mu) + t(1 - y)%*%log(1-mu)
  loglik <- sum(logl)
  #return(-loglik)  # package minimizes objective function so we return neg. log-likelihood
  return(ifelse(is.nan(loglik), Large, -loglik))
}

##### Linear Constraint (No Z)
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
lincon1_z <- function(param, q, y, x.s){
  K <- length(x.s)
  i <- 1
  start <- K+2
  end <- K+1+x.s[i]
  ret <- rep(0, K)
  while(i <= K){
    weights <- param[start:end]
    assign(paste0("sum", i), sum(weights))
    ret[i] <- eval(parse(text=paste0("sum", i)))
    i <- i + 1
    start <- end + 1
    end <- end + x.s[i]
  }
  return(ret)
}

##### Inequality Constraints (No Z)
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
ineq_z <- function(param, q, y, x.s){
  K <- length(x.s)
  i <- 1
  start <- K+2
  end <- K+1+x.s[i]
  ret.b <- rep(0, K)
  ret.w <- numeric()
  while(i <= K){
    assign(paste0("b", i), param[(i+1)])
    assign(paste0("weights", i), param[start:end])
    ret.b[i] <- eval(parse(text=paste0("b", i)))
    ret.w <- c(ret.w, eval(parse(text=paste0("weights", i))))
    i <- i + 1
    start <- end + 1
    end <- end + x.s[i]
  }
  return(c(ret.b, ret.w))
}

##### WQS Fit Internal (No Z)
#' title
#'
#' description
#'
#' deets
#'
#' @noRd
wqs.fit.internal_z <- function(q, y, w, x.s, func){
  K <- length(x.s)
  i <- 1
  start <- 1
  end <- x.s[i]
  #mod <- "b0"
  WQS <- numeric()
  temp.n <- "y"
  while(i <= K){
    assign(paste0("w", i), w[start:end])
    assign(paste0("q", i), q[,start:end])
    assign(paste0("WQS", i), eval(parse(text=paste0("q", i))) %*% eval(parse(text=paste0("w", i))))
    WQS <- cbind(WQS, eval(parse(text=paste0("WQS", i))))
    temp.n <- c(temp.n, paste0("WQS", i))
    #mod <- paste(mod, "+", paste0("b", i),"*", paste0("q", i), "%*%", paste0("w", i), sep="")
    i <- i + 1
    start <- end + 1
    end <- end + x.s[i]
  }
  temp <- data.frame(cbind(y, WQS))
  names(temp) <- temp.n
  if (func == "continuous"){
    fit <- glm2(y ~ ., data = temp, family = "gaussian"(link = identity))
  } else if (func == "binary") {
    fit <- glm2(y ~ ., data = temp, family = "binomial"(link = logit))
  }
  #fit <- lm(y ~ ., data = temp)

  out <- list(temp[,-1], fit)
  names(out) <- c("WQS", "fit")
  return(out)
}

#' ######## Main Function ######## (estimates weights across bootstrap samples for WQS Regression)
#' #' @export
#' wqs.fit <- function(y, y.train = NULL, x, x.train = NULL, z = NULL, z.train = NULL, x.s, B=100, n.quantiles=4, pars = NULL, func,
#'                     eqB = NULL,
#'                     ineqLB = NULL,
#'                     ineqUB = NULL,
#'                     tol = 1e-6,
#'                     delta = 1e-6){
#'
#'   #if (func != "continuous" || "binary")
#'
#'    K <- length(x.s) # Total number of component groups
#'    C <- dim(x)[2] # Total number of components
#'
#' ### Defining pars, function, and all training variables if left as NULL, for Z and no Z cases
#'
#'    if (is.null(z)){         # No Z case
#'                 print("CP no Z setup") # checkpoint
#'            if (is.null(pars)){    # Defining the vector of initial values if pars left as NULL
#'              w_inits <- numeric()
#'              for (i in 1:K){
#'                w_inits <- c(w_inits, rep(1/x.s[i], x.s[i]))
#'              }
#'              pars <- c(rep(0, K+1), w_inits)
#'            }
#'
#'            if (length(pars) != sum(K+1, C)){ stop("The length of pars much equal the total number of betas and weights")}
#'
#'
#'            if (is.null(y.train)){   # Setting all data as training if _.train variables left as NULL
#'              y.train <- y
#'              x.train <- x
#'            }
#'
#'            if (func == "continuous"){     # Setting Objective Function
#'              fun <- cont_z
#'            } else if (func == "binary"){
#'              fun <- bin_z
#'            }
#'
#'    } else {                # Z case
#'           print("CP Z setup")       # checkpoint
#'         p <- dim(z)[2] # Total number of covariates
#'
#'            if (is.null(pars)){        # Defining the vector of initial values if pars left as NULL
#'              w_inits <- numeric()
#'                 for (i in 1:K){
#'                   w_inits <- c(w_inits, rep(1/x.s[i], x.s[i]))
#'                 }
#'              pars <- c(rep(0, K+1), w_inits, rep(.1, p))
#'            }
#'
#'            if (length(pars) != sum(K+1, C, p)){ stop("The length of pars much equal the total number of betas, weights, and covariates")}
#'
#'           if (is.null(y.train)){   # Setting all data as training if _.train variables left as NULL
#'             y.train <- y
#'             x.train <- x
#'             z.train <- z
#'           }
#'
#'           if (func == "continuous"){    # Setting Objective Function
#'             print("CP cont obj fn") # checkpoint
#'             fun <- cont
#'           } else if (func == "binary"){
#'             print("CP binary obj fn") # checkpoint
#'             fun <- bin
#'           }
#'
#'     }
#' ### Defining eqB variable if left as NULL
#'    if (is.null(eqB)){
#'      eqB <- rep(1, K)
#'    }
#'    if (length(eqB) != K){ stop("The length of eqB must equal the number of component groups")}
#' ### Defining ineqLB variable if left as NULL
#'    if (is.null(ineqLB)){
#'      ineqLB <- rep(-2, K)
#'    }
#'    if(length(ineqLB) != K){ stop("The length of ineqLB must equal the number of component groups")}
#'    ineqLB <- c(ineqLB, rep(0,C))
#' ### Defining ineqUB variable if left as NULL
#'    if (is.null(ineqUB)){
#'      ineqUB <- rep(2, K)
#'    }
#'    if(length(ineqUB) != K){ stop("The length of ineqUB must equal the number of component groups")}
#'    ineqUB <- c(ineqUB, rep(1,C))
#' ### Taking quantiles
#'   q <- quantile.fn(x, n.quantiles)
#'   q.train <- quantile.fn(x.train, n.quantiles)
#'
#' ###
#'   if (B > 1) {
#'
#'             if (is.null(z)){
#'                 print("CP No Z, B>1") # checkpoint
#'               result <- matrix(0, nrow = B, ncol = length(pars))
#'               colnames(result) <- names(pars)
#'               convergence <- rep(0, B) #0 indicates convergence; 1 or 2 indicates non-convergence
#'
#'               betas <- numeric()
#'               tests <- numeric()
#'               pvals <- numeric()
#'
#'               #---------------------------- BOOTSTRAP ROUTINE -----------------------------------
#'               for (b in 1:B) {
#'                 # draw random sample (of same size as training data set) with replacement
#'                 samp <- sample(1:length(y.train), replace = TRUE)
#'
#'                 y.b <- as.vector(y.train[samp])
#'                 q.b <- q.train[samp,]
#'
#'                 result.b <- solnp(pars, fun, eqfun = lincon1_z, eqB, ineqfun = ineq_z, ineqLB, ineqUB, LB=NULL, UB=NULL, q.b, y.b, x.s,
#'                                   control = list(tol = tol,delta = delta, trace = 0))
#'
#'                 result[b,] <- result.b$pars
#'                 convergence[b] <- result.b$convergence
#'
#'                 w <- result.b$pars[(K+2):(K+1 + sum(x.s))]
#'                 fit <- wqs.fit.internal_z(q.b, y.b, w, x.s, func)$fit
#'                 b.b <- fit$coefficients[-1]
#'                 t.b <- summary(fit)$coefficients[-1,3]
#'                 p.b <- summary(fit)$coefficients[-1,4]
#'                 betas <- rbind(betas, b.b)
#'                 tests <- rbind(tests, t.b)
#'                 pvals <- rbind(pvals, p.b)
#'               }
#'
#'               out1 <- list(result, convergence, betas, tests, pvals)
#'               print(out1)
#'               names(out1) <- c("params", "convergence", "betas", "tests", "pvals")
#'
#'
#'               wts.matrix <- out1$params[, (K+2):(K+1+ sum(x.s))]
#'               test_stat <- out1$tests
#'               print("CP teststat run") # checkpoint
#'               weights <- teststat.fn2(wts.matrix, test_stat, x.s)
#'
#'               final <- wqs.fit.internal_z(q, y, weights, x.s, func)
#'               print("CP final fit done") # checkpoint
#'               return(final)
#'
#'               } else {
#'
#'               print("CP Z, B>1") #checkpoint
#'               # initialize matrix for parameter estimates (from estimation step)
#'               result <- matrix(0, nrow = B, ncol = length(pars))
#'               colnames(result) <- names(pars)
#'               convergence <- rep(0, B) #0 indicates convergence; 1 or 2 indicates non-convergence
#'
#'               betas <- numeric()
#'               tests <- numeric()
#'               pvals <- numeric()
#'
#'                 #---------------------------- BOOTSTRAP ROUTINE -----------------------------------
#'                 for (b in 1:B) {
#'                   # draw random sample (of same size as training data set) with replacement
#'                   samp <- sample(1:length(y.train), replace = TRUE)
#'
#'                   y.b <- as.vector(y.train[samp])
#'                   q.b <- q.train[samp,]
#'                   z.b <- as.matrix(z.train[samp,])
#'                   rownames(z.b) <- NULL
#'
#'                   result.b <- solnp(pars, fun, eqfun = lincon1, eqB, ineqfun = ineq, ineqLB, ineqUB, q.b, z.b, y.b, x.s, LB = NULL, UB = NULL,
#'                                     control = list(tol = tol, delta = delta, trace = 0))
#'
#'                   result[b,] <- result.b$pars
#'                   convergence[b] <- result.b$convergence
#'
#'                   w <- result.b$pars[(K+2):(K+1 + sum(x.s))] # w: estimated weights from training data (boot strap samples)
#'                   fit <- wqs.fit.internal(q.b, z.b, y.b, w, x.s, func)$fit
#'                   # Remove first and last coefficients for intercept and Z - assume only 1 Z for now
#'                   b.b <- fit$coefficients[-c(1,(K+2))]
#'                   t.b <- summary(fit)$coefficients[-c(1,(K+2)),3]
#'                   p.b <- summary(fit)$coefficients[-c(1,(K+2)),4]
#'                   betas <- rbind(betas, b.b)
#'                   tests <- rbind(tests, t.b)
#'                   pvals <- rbind(pvals, p.b)
#'                 }
#'                 print("CP - bootstrap over") #checkpoint
#'               out1 <- list(result, convergence, betas, tests, pvals)
#'               print(out1)
#'               print("CP out1 created") # checkpoint
#'               names(out1) <- c("params", "convergence", "betas", "tests", "pvals")
#'               #return(out)
#'
#'               wts.matrix <- out1$params[, (K+2):(K+1+ sum(x.s))]
#'               test_stat <- out1$tests
#'               print("CP teststat run?") # checkpoint
#'               weights <- teststat.fn2(wts.matrix, test_stat, x.s)
#'               print("CP teststat run") # checkpoint
#'
#'               final <- wqs.fit.internal(q, z, y, weights, x.s, func)
#'               print("CP last fit run") # checkpoint
#'               return(final)
#'             }
#'
#'   } else if (B == 1) {
#'
#'               if (is.null(z)){
#'
#'                 result <- solnp(pars, fun, eqfun = lincon1_z, eqB,
#'                                 ineqfun = ineq_z, ineqLB, ineqUB, LB = NULL, UB = NULL, q = q.train, y = y.train, x.s = x.s,
#'                                 control = list(tol = tol, delta = delta, trace = 0))
#'
#'                 par.estimates <- result$pars
#'
#'                 w_names <- character()
#'                 for (i in 1:K){
#'                   w_names <- c(w_names, paste0("w", i, 1:x.s[i]))
#'                 }
#'                 #names(par.estimates) <- c("b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", paste0("w1", 1:x.s[1]), paste0("w2", 1:x.s[2]), paste0("w3", 1:x.s[3]), paste0("w4", 1:x.s[4]), paste0("w5", 1:x.s[5]), paste0("w6", 1:x.s[6]), paste0("w7", 1:x.s[7]), paste0("z", 1:p))
#'                 names(par.estimates) <- c(paste0("b", 0:K), w_names)
#'                 round(par.estimates, 2)
#'                 weights <- par.estimates[(K+2):(K+1+sum(x.s))]
#'
#'                 final <- wqs.fit.internal_z(q, y, weights, x.s, func)
#'                 return(final)
#'
#'               } else {
#'
#'                 result <- solnp(pars, fun, eqfun = lincon1, eqB,
#'                                 ineqfun = ineq, ineqLB, ineqUB, LB = NULL, UB = NULL, q = q.train, y = y.train, x.s = x.s,
#'                                 z = z.train, control = list(tol = tol, delta = delta, trace = 0))
#'
#'                 par.estimates <- result$pars
#'
#'                 w_names <- character()
#'                 for (i in 1:K){
#'                   w_names <- c(w_names, paste0("w", i, 1:x.s[i]))
#'                 }
#'                 #names(par.estimates) <- c("b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", paste0("w1", 1:x.s[1]), paste0("w2", 1:x.s[2]), paste0("w3", 1:x.s[3]), paste0("w4", 1:x.s[4]), paste0("w5", 1:x.s[5]), paste0("w6", 1:x.s[6]), paste0("w7", 1:x.s[7]), paste0("z", 1:p))
#'                 names(par.estimates) <- c(paste0("b", 0:K), w_names, paste0("z", 1:p))
#'                 round(par.estimates, 2)
#'                 weights <- par.estimates[(K+2):(K+1+sum(x.s))]
#'
#'                 final <- wqs.fit.internal(q, z, y, weights, x.s, func)
#'                 return(final)
#'               }
#'     } else {
#'
#'       stop("B parameter must have positive, nonzero value")
#'
#'     }
#' }
