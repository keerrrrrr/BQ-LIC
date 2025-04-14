# Bayesian Quantile Regression Model with Linear Inequality Constraints

This code is designed to solve Bayesian quantile regression model with linear inequality constraints. Leveraging asymmetric Laplace distributions, we propose two novel Gibbs sampling methods based on truncated normal and truncated Laplace prior distributions, respectively.

## Packages
### Install the `R` package `BQ-LIC` from [GitHub](https://github.com/keerrrrrr/BQ-LIC)
```R
# install.packages('devtools')
devtools::install_github("keerrrrrr/BQ-LIC")
```

## Example
This is a basic example which shows you how to use the functions in this package to solve the Bayesian quantile regression model with linear inequality constraints.

```R
library(GIGrvg)
library(MASS)
library(invgamma)
library(tmvmixnorm)
library(mvtnorm)
library(BQ-LIC)


# Number of observations
n <- 200

# Number of predictors
d <- 3

#parameter
beta <- c(0.3,0.7,-0.5)

#constraints
R <- matrix(c(1,0,0,1,0,1),2,3)
b <- c(0.2,0)

#quantile
p <- 0.5

#data
x <- rmvnorm(n,rep(0,d),diag(d))
epsilon <- rnorm(n,0,1)
epsilon1 <- epsilon-qnorm(p)
y <- x%*%beta+epsilon1

#initial value
intsp <- c(0.5,0.5,-0.1)


#parameter estimation for the Bayesian quantile regression model with linear inequality constraints based on the truncated normal prior
beta_LIC_TN <- BQ_LIC_TN(x, y, p, R, b)

#parameter estimation for the Bayesian quantile regression model with linear inequality constraints based on the truncated Laplace prior
beta_LIC_TL <- BQ_LIC_TL(x, y, p, R, b)
```
