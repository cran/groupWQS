#' Grouped WQS Regression
#'
#' This function fits a grouped weighted quantile sum (GWQS) regression model.
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
#' @param pars A vector of initial values, listed in order: beta naught intercept and group index beta coefficients, individual chemical weight coefficients, and covariate coefficients.
#' @param func The objective function to be used (must match outcome data type); currently only fun args "continuous" or "binary" are supported.
#' @param ineqLB Vector of lower bounds for betas and weights, set to -2 by default.
#' @param ineqUB Vector of upper bounds for betas and weights, set to 2 be default.
#' @param tol Tolerance level for bootstrap convergence.
#' @param delta Step size for bootstrap procedure.
#'
#' @return A list of 3 containing the GWQS estimate based on calculated weights, the GWQS model fit to validation data, and weight estimates
#'
#' @examples
#' data("WQSdata")
#' group_list <- list(c("X1", "X2", "X3"), c("X4", "X7"), c("X5", "X6", "X9", "X8"))
#' x.s <- make.x.s(WQSdata, 3, group_list)
#' X <- make.X(WQSdata, 3, group_list)
#' Y <- WQSdata$y
#' results <- gwqs.fit(y = Y, x = X, x.s = x.s, B=1, func = "continuous")
#'
#' @import stats
#' @import Rsolnp
#' @import glm2
#'
#' @export
gwqs.fit <- function(y, y.train = NULL, x, x.train = NULL, z = NULL, z.train = NULL, x.s, B=100, n.quantiles=4, pars = NULL, func,
                    ineqLB = NULL,
                    ineqUB = NULL,
                    tol = 1e-6,
                    delta = 1e-6){

 #    if (any(is.NA(y.train))){
 #      stop("NAs in y.train")
 #    }
 #    if (any(is.NA(x.train))){
 #        stop("NAs in x.train")
 #    }
 #
 #  if (!is.null(y.train)){
 #      # Check training data
 #      check_train <- check_xyz(xdata=x.train, ydata=y.train, zdata=z.train)
 #      if(check_train[1]) {
 #          stop("x.train must be a matrix")
 #      }
 #      if(check_train[2]) {
 #          stop("check dimensions of training data2")
 #      }
 #      if(check_train[3]) {
 #          stop("check dimensions of training data3")
 #      }
 # }

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

      z <- as.matrix(z)
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
  # Renaming weight values, turning the weight vector to a list, adding to returned object "final"
      names(weights) <- colnames(x)
      wlist <- list()
      count <- 1
      for (i in 1:length(x.s)){
        wlist[[i]] <- weights[count:(count+x.s[i]-1)]
        count <- count + x.s[i]
      }
      final[[3]] <- wlist

      names(final) <- c("GWQS", "fit", "weights")
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
      # Renaming weight values, turning the weight vector to a list, adding to returned object "final"
      names(weights) <- colnames(x)
      wlist <- list()
      count <- 1
      for (i in 1:length(x.s)){
        wlist[[i]] <- weights[count:(count+x.s[i]-1)]
        count <- count + x.s[i]
      }
      final[[3]] <- wlist

      names(final) <- c("GWQS", "fit", "weights")
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
      # Renaming weight values, turning the weight vector to a list, adding to returned object "final"
      names(weights) <- colnames(x)
      wlist <- list()
      count <- 1
      for (i in 1:length(x.s)){
        wlist[[i]] <- weights[count:(count+x.s[i]-1)]
        count <- count + x.s[i]
      }
      final[[3]] <- wlist

      names(final) <- c("GWQS", "fit", "weights")
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
      # Renaming weight values, turning the weight vector to a list, adding to returned object "final"
      names(weights) <- colnames(x)
      wlist <- list()
      count <- 1
      for (i in 1:length(x.s)){
        wlist[[i]] <- weights[count:(count+x.s[i]-1)]
        count <- count + x.s[i]
      }
      final[[3]] <- wlist

      names(final) <- c("GWQS", "fit", "weights")
      return(final)
    }
  } else {

    stop("B parameter must have positive, nonzero value")

  }
}

