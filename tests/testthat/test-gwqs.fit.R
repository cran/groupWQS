library(testthat)
library(groupWQS)

context("gwqs.fit() tests")

# Prepping Covariate Data

set.seed(1990)
female <- rbinom(1000, 1,.5)
age <- abs(rnorm(1000, mean=0, sd=5))
Z <- cbind(female, age)
Z_bin <- Z[,1]

# Prepping Data with Binary Outcome

data("simdata")

group_list1 <- list(c("pcb_118", "pcb_138", "pcb_153"),
                   c("as", "cu"),
                   c("carbaryl", "propoxur", "methoxychlor", "diazinon", "chlorpyrifos"),
                   c("sn"))
x.s1 <- make.x.s(simdata, 4, group_list1)
X1 <- make.X(simdata, 4, group_list1)
Y1 <- simdata$Y

# Splitting Binary Data

smp_size <- floor(0.50 * nrow(X1))
set.seed(1990)
train_ind <- sample(seq_len(nrow(X1)), size = smp_size)

X1.train <- X1[train_ind, ]
X1.valid <- X1[-train_ind, ]

Z.train <- Z[train_ind,]
Z.valid <- Z[-train_ind,]

Y1.train <- Y1[train_ind]
Y1.valid <- Y1[-train_ind]

# Prepping Data with Continuous Outcome

data("WQSdata")

group_list2 <- list(c("X1", "X2", "X3"),
                    c("X4"),
                    c("X5","X6"))
x.s2 <- make.x.s(WQSdata, 3, group_list2)
X2 <- make.X(WQSdata, 3, group_list2)
Y2 <- WQSdata$y


# Splitting Continuous Data

smp_size2 <- floor(0.50 * nrow(X2))
set.seed(1990)
train_ind2 <- sample(seq_len(nrow(X2)), size = smp_size2)

X2.train <- X2[train_ind2, ]
X2.valid <- X2[-train_ind2, ]

Z2.train <- Z[train_ind2,]
Z2.valid <- Z[-train_ind2,]

Y2.train <- Y2[train_ind2]
Y2.valid <- Y2[-train_ind2]

# Binary Tests

test_that("(1) binary, no covariates, no split data, B>1 ",{

    expect_error(gwqs.fit(y = Y1, x = X1, x.s = x.s1, B=5, n.quantiles = 4, func = "binary"), NA)
})

test_that("(2) binary, no covariates, no split data, B=1 ",{
    expect_error(gwqs.fit(y = Y1, x = X1, x.s = x.s1, B=1, n.quantiles = 4, func = "binary"), NA)
})

test_that("(3) binary, covariates, no split data, B>1 ",{

    expect_error(gwqs.fit(y = Y1, x = X1, z = Z, x.s = x.s1, B=5, n.quantiles = 4, func = "binary"), NA)
})

test_that("(4) binary, single covariate, no split data, B>1 ",{

    expect_error(gwqs.fit(y = Y1, x = X1, z = Z_bin, x.s = x.s1, B=5, n.quantiles = 4, func = "binary"), NA)
})

test_that("(5) binary, covariates, no split data, B=1 ",{
    expect_error(gwqs.fit(y = Y1, x = X1, z = Z, x.s = x.s1, B=1, n.quantiles = 4, func = "binary"), NA)
})

test_that("(6) binary, single covariates, no split data, B=1 ",{
    expect_error(gwqs.fit(y = Y1, x = X1, z = Z_bin, x.s = x.s1, B=1, n.quantiles = 4, func = "binary"), NA)
})


test_that("(7) binary, covariates, split data, B>1",{

    expect_error(gwqs.fit(y = Y1.valid, y.train = Y1.train, x = X1.valid, x.train = X1.valid, z = Z.valid, z.train = Z.train,
                          x.s = x.s1, B=5, n.quantiles = 4, func = "binary"), NA)
})

test_that("(8) binary, covariates, split data, B=1",{

    expect_error(gwqs.fit(y = Y1.valid, y.train = Y1.train, x = X1.valid, x.train = X1.valid, z = Z.valid, z.train = Z.train,
                          x.s = x.s1, B=1, n.quantiles = 4, func = "binary"), NA)
})

test_that("(9) binary, no covariates, split data, B>1",{

    expect_error(gwqs.fit(y = Y1.valid, y.train = Y1.train, x = X1.valid, x.train = X1.valid,
                          x.s = x.s1, B=5, n.quantiles = 4, func = "binary"), NA)
})

test_that("(10) binary, no covariates, split data, B=1",{

    expect_error(gwqs.fit(y = Y1.valid, y.train = Y1.train, x = X1.valid, x.train = X1.valid,
                          x.s = x.s1, B=1, n.quantiles = 4, func = "binary"), NA)
})

# Continuous Tests

test_that("(11) continuous, no covariates, no split data, B>1 ",{

    expect_error(gwqs.fit(y = Y2, x = X2, x.s = x.s2, B=5, n.quantiles = 4, func = "continuous"), NA)
})

test_that("(12) continuous, no covariates, no split data, B=1 ",{
    expect_error(gwqs.fit(y = Y2, x = X2, x.s = x.s2, B=1, n.quantiles = 4, func = "continuous"), NA)
})

test_that("(13) continuous, covariates, no split data, B>1 ",{

    expect_error(gwqs.fit(y = Y2, x = X2, z = Z, x.s = x.s2, B=5, n.quantiles = 4, func = "continuous"), NA)
})

test_that("(14) continuous, single covariate, no split data, B>1 ",{

    expect_error(gwqs.fit(y = Y2, x = X2, z = Z_bin, x.s = x.s2, B=5, n.quantiles = 4, func = "continuous"), NA)
})

test_that("(15) continuous, covariates, no split data, B=1 ",{
    expect_error(gwqs.fit(y = Y2, x = X2, z = Z, x.s = x.s2, B=1, n.quantiles = 4, func = "continuous"), NA)
})

test_that("(16) continuous, single covariates, no split data, B=1 ",{
    expect_error(gwqs.fit(y = Y2, x = X2, z = Z_bin, x.s = x.s2, B=1, n.quantiles = 4, func = "continuous"), NA)
})


test_that("(17) continuous, covariates, split data, B>1",{

    expect_error(gwqs.fit(y = Y2.valid, y.train = Y2.train, x = X2.valid, x.train = X2.valid, z = Z2.valid, z.train = Z2.train,
                          x.s = x.s2, B=5, n.quantiles = 4, func = "continuous"), NA)
})

test_that("(18) continuous, covariates, split data, B=1",{

    expect_error(gwqs.fit(y = Y2.valid, y.train = Y2.train, x = X2.valid, x.train = X2.valid, z = Z2.valid, z.train = Z2.train,
                          x.s = x.s2, B=1, n.quantiles = 4, func = "continuous"), NA)
})

test_that("(19) continuous, no covariates, split data, B>1",{

    expect_error(gwqs.fit(y = Y2.valid, y.train = Y2.train, x = X2.valid, x.train = X2.valid,
                          x.s = x.s2, B=5, n.quantiles = 4, func = "continuous"), NA)
})

test_that("(20) continuous, no covariates, split data, B=1",{

    expect_error(gwqs.fit(y = Y2.valid, y.train = Y2.train, x = X2.valid, x.train = X2.valid,
                          x.s = x.s2, B=1, n.quantiles = 4, func = "continuous"), NA)
})
