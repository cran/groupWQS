#' Generates Plots of weights by group
#'
#' This function takes the object created by the wqs.fit function and a vector of group names and generates a
#' random forest variable importance plot for each group. The weights in each group are listed in descending order.
#'
#' @param fit.object The object that is returned by the wqs.fit function
#' @param group.names A string vector containing the name of each group included in the GWQS regression. Will be used for plot titles.
#'
#' @return A plot for each group of the GWQS regression
#'
#' @examples
#' data("WQSdata")
#' group_list <- list(c("X1", "X2", "X3"), c("X4", "X7"), c("X5", "X6", "X9", "X8"))
#' chem_groups <- c("PCBs", "Metals", "Insecticides")
#' x.s <- make.x.s(WQSdata, 3, group_list)
#' X <- make.X(WQSdata, 3, group_list)
#' Y <- WQSdata$y
#' results <- gwqs.fit(y = Y, x = X, x.s = x.s, B=1, func = "continuous")
#' weight.plot(results, chem_groups)
#' @export
weight.plot <- function(fit.object, group.names){
    plot_num <- length(fit.object[[3]])

    for (i in 1:plot_num){
        weight_vec <- fit.object[[3]][[i]]
        sorted_vec <- sort(weight_vec)
        graphics::dotchart(sorted_vec, labels = names(sorted_vec), xlab = "Weight", pch = 19)
        graphics::title(main = paste("Weights for", group.names[i]), cex.main = 1.00)
    }
}



