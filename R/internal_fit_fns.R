### With Covariates (Z) ###

wqs.fit.internal <- function(q, z, y, w, x.s, func){
    K <- length(x.s)
    if (is.vector(z)){
        p <- 1
    } else {
        p <- dim(z)[2] # Total number of covariates
    }
    i <- 1
    start <- 1
    end <- x.s[i]
    #mod <- "b0"
    GWQS <- numeric()
    temp.n <- "y"
    while(i <= K){
        assign(paste0("w", i), as.matrix(w[start:end]))
        assign(paste0("q", i), q[,start:end])
        assign(paste0("GWQS", i), eval(parse(text=paste0("q", i))) %*% eval(parse(text=paste0("w", i))))
        GWQS <- cbind(GWQS, eval(parse(text=paste0("GWQS", i))))
        temp.n <- c(temp.n, paste0("GWQS", i))
        #mod <- paste(mod, "+", paste0("b", i),"*", paste0("q", i), "%*%", paste0("w", i), sep="")
        i <- i + 1
        start <- end + 1
        end <- end + x.s[i]
    }
    temp <- data.frame(cbind(y, GWQS, z))
    temp.n <- c(temp.n, paste0("Z", seq(1, p, 1)))
    names(temp) <- temp.n
    if (func == "continuous"){
        fit <- glm2(y ~ ., data = temp, family = "gaussian"(link = identity))
    } else if (func == "binary") {
        fit <- glm2(y ~ ., data = temp, family = "binomial"(link = logit))
    }
    #fit <- lm(y ~ ., data = temp)
    
    out <- list(temp[,-1], fit)
    names(out) <- c("GWQS", "fit")
    return(out)
}

### Without Covariates (no Z) ###

wqs.fit.internal_z <- function(q, y, w, x.s, func){
    K <- length(x.s)
    i <- 1
    start <- 1
    end <- x.s[i]
    #mod <- "b0"
    GWQS <- numeric()
    temp.n <- "y"
    while(i <= K){
        assign(paste0("w", i), as.matrix(w[start:end]))
        assign(paste0("q", i), q[,start:end])
        assign(paste0("GWQS", i), eval(parse(text=paste0("q", i))) %*% eval(parse(text=paste0("w", i))))
        GWQS <- cbind(GWQS, eval(parse(text=paste0("GWQS", i))))
        temp.n <- c(temp.n, paste0("GWQS", i))
        #mod <- paste(mod, "+", paste0("b", i),"*", paste0("q", i), "%*%", paste0("w", i), sep="")
        i <- i + 1
        start <- end + 1
        end <- end + x.s[i]
    }
    temp <- data.frame(cbind(y, GWQS))
    names(temp) <- temp.n
    if (func == "continuous"){
        fit <- glm2(y ~ ., data = temp, family = "gaussian"(link = identity))
    } else if (func == "binary") {
        fit <- glm2(y ~ ., data = temp, family = "binomial"(link = logit))
    }
    #fit <- lm(y ~ ., data = temp)
    
    out <- list(temp[,-1], fit)
    names(out) <- c("GWQS", "fit")
    return(out)
}