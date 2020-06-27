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