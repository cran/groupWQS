### With Covariates (Z) ###

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

### Without Covariates (no Z) ###

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