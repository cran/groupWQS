### With Covariates (Z) ###

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
        assign(paste0("w", i), as.matrix(param[w.start:w.end]))
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
    #theta <- param[C]   # assume only one covariate for now
    theta <- as.matrix(param[(sum(x.s)+K+2):(length(param))])
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

### Without Covariates (no Z) ###

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
        assign(paste0("w", i), as.matrix(param[w.start:w.end]))
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