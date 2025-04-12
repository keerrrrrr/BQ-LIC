#' @title Bayesian inference of quantile regression models with linear inequality constraints based on truncated Laplace prior
#'
#' @param x data X, n*d matrix
#' @param y data Y, n*1 vector
#' @param p quantile level
#' @param R constraint matrix with the dimension m*d
#' @param b constraint vector with the dimension m
#'
#' @return estimated parameter, 1*d vector
#' @export
#'
#' @examples
#' n <- 200
#' d <- 3
#' beta <- c(0.3,0.7,-0.5)
#' R <- matrix(c(1,0,0,1,0,1),2,3)
#' b <- c(0.2,0)
#' p <- 0.5
#' x <- rmvnorm(n,rep(0,d),diag(d))
#' epsilon <- rnorm(n,0,1)
#' epsilon1 <- epsilon-qnorm(p)
#' y <- x%*%beta+epsilon1
#' intsp <- c(0.5,0.5,-0.1) #initial value
#' beta_LIC_TN <- BQ_LIC_TN(x,y,p,R,b)

BQ_LIC_TN <- function(x,y,p,R,b){
  iter_total <- 5000
  n <- length(y)
  d <- ncol(x)
  m <- length(b)
  mu0 <- rep(0,d)
  Sigma0 <- 100*diag(d)
  sigma <- 1
  n0 <- 3
  s0 <- 0.01
  beta_samples <- matrix(0,3000,d)
  theta <- (1-2*p)/(p*(1-p))
  tau2 <- 2/(p*(1-p))
  beta_p <- rep(1,d)
  z <- rexp(n,1)
  beta_p_record <- rep(list(), iter_total)
  iter <- 1
  for (iter in 1:iter_total){
    ## Sigma_hat
    Sigma_hat <- solve(Sigma0)
    for (i in 1:n){
      Sigma_hat<- Sigma_hat + (x[i,]%*%t(x[i,]) / (tau2*z[i]*sigma))
    }
    Sigma_hat <- solve(Sigma_hat)

    ## mu_hat
    mu_hat <- solve(Sigma0)%*%mu0
    for (i in 1:n){
      mu_hat <- mu_hat + matrix(x[i,],d,1)*(y[i]-theta*z[i]) / (tau2*z[i]*sigma)
    }
    mu_hat <- Sigma_hat %*% mu_hat
    beta_p <- rtmvn(n=1, Mean=mu_hat, Sigma_hat, D=R, lower=b,upper=rep(Inf,m),int=intsp, burn=10)
    beta_p_record[[iter]] <- beta_p

    ## z
    for (i in 1:n){
      delta_i_hat <- (y[i]-(t(x[i,])%*%beta_p))^2/(tau2*sigma)
      gamma_i_hat <- (theta^2/tau2 +2)/sigma
      z[i] <- rgig(1, 0.5,delta_i_hat,gamma_i_hat)
    }

    ##sigma
    n_tidle <- 2*n0+3*n
    s_tidle <- 2*s0+2*sum(z)+sum((y-x%*%beta_p-theta*z)^2/(tau2*z))
    sigma <- rinvgamma(1,n_tidle/2,s_tidle/2)
  }
  for(i in 1:d){
    beta_samples[,i] <- unlist(lapply(beta_p_record, function(x){return(x[i])}))[2001:iter_total]
  }
  return(colMeans(beta_samples))
}

