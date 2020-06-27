#' Forms matrix of components
#'
#' This function returns a matrix of component variables, X. The user can specify
#' the desired chemicals and order by creating a list of string vectors, each vector containing the variable names of all
#' desired elements of that group.
#'
#'
#' @param df A dataframe containing named component variables
#' @param num.groups An integer representing the number of component groups desired
#' @param groups A list, each item in the list being a string vector of variable names for one component group
#'
#' @return A matrix of component variables
#'
#' @examples
#' data("WQSdata")
#' group_list <- list(c("X1", "X2", "X3"), c("X4", "X7"), c("X5", "X6", "X9", "X8"))
#' X <- make.X(WQSdata, 3, group_list)
#' X
#' @export
make.X <- function(df, num.groups, groups){
    if (any(class(df)=="spec_tbl_df")){df <- as.data.frame(df)}
 xrow <- dim(df)[1]
 x.s <- numeric()
   for (i in 1:num.groups){
          x.s[i] <- length(groups[[i]])
        }
 xcol <- sum(x.s)

 X <- matrix(0, nrow=xrow, ncol=xcol)
 counter <-1

  for (j in 1:num.groups){
        for (k in 1:x.s[j]){
          x <- groups[[j]][k]
          X[, counter] <- df[, x]
          counter <- counter + 1
        }
    counter <- counter
  }
 colnames(X) <- unlist(groups)
return(X)
}


