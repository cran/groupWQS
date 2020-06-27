### With Covariates (Z) ###

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
        assign(paste0("w", i), as.matrix(param[w.start:w.end]))
        assign(paste0("q", i), q[,q.start:q.end])
        mod <- paste(mod, "+", paste0("b", i),"*", paste0("q", i), "%*%", paste0("w", i), sep="")
        i <- i + 1
        w.start <- w.end + 1
        w.end <- w.end + x.s[i]
        q.start <- q.end + 1
        q.end <- q.end + x.s[i]
    }
    
    # Add covariate to model
    
    theta <- as.matrix(param[(sum(x.s)+K+2):(length(param))])
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

### Without Covariates (no Z) ###

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
        assign(paste0("w", i), as.matrix(param[w.start:w.end]))
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