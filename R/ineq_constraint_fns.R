### With Covariates (Z) ###

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

### Without Covariates ###

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