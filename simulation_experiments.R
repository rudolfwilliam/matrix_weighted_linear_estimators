######## Obtain an estimate of the expected risk of different weighting methods

source("./data_generation.R")
source("./weighting_related_work.R")
library(glmnet)
library(MASS)

set.seed(1)

############ SETUP CODE (you may change this parameters as you wish) ############

# whether the entries of X should be dependent on each other
dependencies <- FALSE

# whether confounding should affect only a few parameters or all of them
sparse <- FALSE

# num iterations for obtaining estimate of expected risk
num.it <- 1000

# initialize SCM parameters
# num features
p <- 30

if(sparse) {
  # num variables that are affected by confounding
  q <- 5
}

# num observational data
n <- 600
# num interventional data
m <- 300

# standard deviations
sd.nz <- 1
sd.nx <- 1
sd.ny <- 1

if(sparse) {
  # only first q features are affected by confounding
  B <- matrix(c(rnorm(q, 0, 1), rep(0, p - q)), ncol=1)
} else {
  B <- matrix(rnorm(p, 0, 1), ncol=1)
}

alpha <- matrix(rnorm(n=p, mean=0, sd=3), ncol=1)
gamma <- 5
# ensure that MSE converges to 0 in the limit
eps <- 1e-04

#################################################################################

if(dependencies) {
  # sample covariance matrix
  mat <- matrix(runif(p^2, min=0, max=sqrt(.5)), p, p)
  # make symmetrical
  cov.mat.x.o <- (mat + t(mat))/2
  diag(cov.mat.x.o) <- 1
  
  # sample covariance matrix
  mat <- matrix(runif(p^2, min=0, max=sqrt(.5)), p, p)
  # make symmetrical
  cov.mat.x.i <- (mat + t(mat))/2
  diag(cov.mat.x.i) <- 1
} else {
  cov.mat.x.o <- diag(rep(sd.nx^2, p))
  cov.mat.x.i <- diag(rep(sd.nx^2, p))
}

# compute sqrt of cov.mat from SVD (required for sampling data)
svd.cov.mat.o <- svd(cov.mat.x.o)
cov.mat.x.o.sqrt <- svd.cov.mat.o$u %*% sqrt(diag(svd.cov.mat.o$d)) %*% t(svd.cov.mat.o$u)

svd.cov.mat.i <- svd(cov.mat.x.i)
cov.mat.x.i.sqrt <- svd.cov.mat.i$u %*% sqrt(diag(svd.cov.mat.i$d)) %*% t(svd.cov.mat.i$u)

# theoretically optimal quantities (unknown in practice)
Delta <- sd.nz^2 * solve(cov.mat.x.o + sd.nz^2 * B %*% t(B)) %*% B * gamma
Cov.X.o <- sd.nz^2 * (B %*% t(B)) + cov.mat.x.o
Cov.Y.X <- B %*% gamma * sd.nz^2 + Cov.X.o %*% alpha
Var.Y <- sd.ny^2 + t(alpha) %*% cov.mat.x.o %*% alpha + sd.nz^2 * t(gamma + t(B) %*% alpha) %*% (gamma + t(B) %*% alpha)
sd.O <- sqrt(Var.Y - t(Cov.Y.X) %*% solve(Cov.X.o) %*% Cov.Y.X)
sd.I <- sqrt(sd.nz^2 * t(gamma) %*% gamma + sd.ny^2)

# initialize error vectors
errs.l2 <- errs.l1 <- errs.standard <- errs.rosenman <- errs.optimal <- errs.I <- errs.O <- errs.pool <- numeric(num.it)

for (it in 1:num.it) {
  # generate data from structural causal model
  data.O <- generate_O(n=n, p=p, sd.nz=sd.nz, sqrt.cov.x=cov.mat.x.o.sqrt, alpha=alpha, gamma=gamma, B=B)
  X.O <- data.O[[1]]
  Y.O <- data.O[[2]]
    
  data.I <- generate_I(n=m, p=p, sd.nz=sd.nz, sqrt.cov.x=cov.mat.x.i.sqrt, alpha=alpha, gamma=gamma, B=B)
  X.I <- data.I[[1]]
  Y.I <- data.I[[2]]
  
  # ground truth covariance matrices (needed for optimal matrix weight)
  Cov.O <- solve(t(X.O) %*% X.O) * as.numeric(sd.O^2)
  Cov.I <- solve(t(X.I) %*% X.I) * as.numeric(sd.I^2)
  
  # OLS estimates and their estimated covariance matrices
  alpha.hat.O.model <- lm(Y.O ~ X.O + 0)
  alpha.hat.I.model <- lm(Y.I ~ X.I + 0)
  alpha.hat.pool.model <- lm(rbind(Y.O, Y.I) ~ rbind(X.O, X.I) + 0)
  
  Cov.O.hat <- solve(t(X.O) %*% X.O) * summary(alpha.hat.O.model)$sigma^2
  Cov.I.hat <- solve(t(X.I) %*% X.I) * summary(alpha.hat.I.model)$sigma^2
  
  # compute bias using ridge regression
  target <- Y.I - X.I %*% alpha.hat.O.model$coefficients
  cv <- cv.glmnet(y=target, x=as.matrix(X.I), alpha=0)
  lambda.min <- cv$lambda.min
  bias.l2 <- -glmnet(y=target, x=X.I, alpha=0, lambda=lambda.min)$beta
  
  # compute bias using Lasso regression
  cv <- cv.glmnet(y=target, x=as.matrix(X.I), alpha=1)
  lambda.min <- cv$lambda.min
  bias.l1 <- glmnet(y=target, x=X.I, alpha=1, lambda=lambda.min)$beta
  
  # bias without ridge regression or Lasso regression
  bias <- alpha.hat.O.model$coefficients - alpha.hat.I.model$coefficients
  
  # weights and resulting models
  weight.l2 <- (Cov.O.hat + bias.l2 %*% t(bias.l2)) %*% solve(Cov.O.hat + Cov.I.hat + bias.l2 %*% t(bias.l2))
  alpha.hat.l2 <- weight.l2 %*% alpha.hat.I.model$coefficients + (diag(1, nrow=p) - weight.l2) %*% alpha.hat.O.model$coefficients
  weight.l1 <- (Cov.O.hat + bias.l1 %*% t(bias.l1)) %*% solve(Cov.O.hat + Cov.I.hat + bias.l1 %*% t(bias.l1))
  alpha.hat.l1 <- weight.l1 %*% alpha.hat.I.model$coefficients + (diag(1, nrow=p) - weight.l1) %*% alpha.hat.O.model$coefficients
  weight.standard <- (Cov.O.hat + bias %*% t(bias)) %*% solve(Cov.O.hat + Cov.I.hat + bias %*% t(bias))
  alpha.hat.standard <- weight.standard %*% alpha.hat.I.model$coefficients + (diag(1, nrow=p) - weight.standard) %*% alpha.hat.O.model$coefficients
  weight.rosenman <- rosenman_weight(alpha.hat.O = alpha.hat.O.model, alpha.hat.I = alpha.hat.I.model, cov.est=Cov.I.hat)
  alpha.hat.rosenman <- weight.rosenman[1] * alpha.hat.I.model$coefficients + (1 - weight.rosenman[1]) * alpha.hat.O.model$coefficients
  weight.optimal <- (Cov.O + Delta %*% t(Delta)) %*% solve(Cov.O + Cov.I + Delta %*% t(Delta))
  alpha.hat.optimal <- weight.optimal %*% alpha.hat.I.model$coefficients + (diag(1, nrow=p) - weight.optimal) %*% alpha.hat.O.model$coefficients
  # errors
  errs.pool[it] <- sum((alpha.hat.pool.model$coefficients - alpha)^2)
  errs.I[it] <- sum((alpha.hat.I.model$coefficients - alpha)^2)
  errs.standard[it] <- sum((alpha.hat.standard - alpha)^2)
  errs.l1[it] <- sum((alpha.hat.l1 - alpha)^2)
  errs.l2[it] <- sum((alpha.hat.l2 - alpha)^2)
  errs.rosenman[it] <- sum((alpha.hat.rosenman - alpha)^2)
  errs.optimal[it] <- sum((alpha.hat.optimal - alpha)^2)
}

print("weighting method, mean squared error, sd of squared error")
print(paste("l2, ", mean(errs.l2), ", ", sd(errs.l2)))
print(paste("l1, ", mean(errs.l1), ", ", sd(errs.l1)))
print(paste("standard/plug-in, ", mean(errs.standard), ", ", sd(errs.standard)))
print(paste("pure interventional, ", mean(errs.I), ", ", sd(errs.I)))
print(paste("pooling, ", mean(errs.pool), ", ", sd(errs.pool)))
print(paste("rosenman, ", mean(errs.rosenman), ", ", sd(errs.rosenman)))
print(paste("optimal (oracle), ", mean(errs.optimal), ", ", sd(errs.optimal)))