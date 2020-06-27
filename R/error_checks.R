check_xyz <- function(xdata,ydata,zdata){
    class_x <- ifelse(class(xdata)== "matrix" | class(xdata)== "data.frame", 0, 1)
    
    n <- length(ydata)
    n.x <- dim(xdata)[1]
    dim_x <- ifelse(n == n.x, 0, 1)
    
    if(is.null(zdata)){dim_z <- 0}
    else{
        if(class(zdata) == "matrix") n.z <- dim(zdata)[1]
        if(class(zdata) == "vector") n.z <- length(zdata)[1]
        dim_z <- ifelse(n != n.z, 1, 0)
    }
    
    return(c(class_x, dim_x, dim_z))
}