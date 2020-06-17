## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(groupWQS)

## ---- include=T---------------------------------------------------------------
data(simdata)
head(simdata)

## ---- include=T---------------------------------------------------------------
group_list <- list(c("pcb_118", "pcb_138",  "pcb_153", "pcb_180", "pcb_192"),
                   c("as", "cu", "pb", "sn"), 
                   c("carbaryl", "propoxur", "methoxychlor", "diazinon", "chlorpyrifos"))
X <- make.X(simdata, num.groups = 3, groups = group_list)
head(X)

## ---- include=T---------------------------------------------------------------
x.s <- make.x.s(simdata, num.groups = 3, groups = group_list)
x.s

## ---- include=T---------------------------------------------------------------
Y <- simdata$Y

## ---- include=T---------------------------------------------------------------
Y.train <- Y[1:700]
Y.valid <- Y[701:1000]
X.train <- X[1:700,]
X.valid <- X[701:1000,]

## ---- include=T---------------------------------------------------------------
results <- gwqs.fit(y = Y.valid, y.train = Y.train, x = X.valid, x.train = X.train, x.s = x.s, 
                   B=100, n.quantiles = 5, func = "binary")

## ---- include=T---------------------------------------------------------------
summary(results$fit)$coefficients

## ---- include=T---------------------------------------------------------------
summary(results$fit)$aic

## ---- include=T---------------------------------------------------------------
exp(confint(results$fit))

## ---- include=T, out.height="60%", out.width="60%", fig.cap=c("Figure 1: Importance plot for chemical concentration weights of GWQS1 index", "Figure 2: Importance plot for chemical concentration weights of GWQS2 index", "Figure 3: Importance plot for chemical concentration weights of GWQS3 index")----
chem_groups <- c("PCBs", "Metals", "Insecticides")
weight.plot(results, chem_groups)

