# helper script to generate confounded and unconfounded data

# confounded data
generate_O <- function(n, p, sd.nz, sqrt.cov.x, alpha, gamma, B) {
  Z <- matrix(rnorm(n=n, sd=sd.nz), ncol=1)
  X <- matrix(rnorm(n=n*p, sd=1, mean=0), ncol=p) %*% sqrt.cov.x + Z %*% t(B)
  Y <- X %*% alpha + gamma * Z + rnorm(n=n)
  
  return(list(X, Y))
}

# unconfounded data
generate_I <- function(n, p, sd.nz, sqrt.cov.x, alpha, gamma, B) {
  Z <- matrix(rnorm(n=n, sd=sd.nz), ncol=1)
  X <- matrix(rnorm(n=n*p, sd=1, mean=0), ncol=p) %*% sqrt.cov.x + Z %*% t(B)
  Z <- matrix(rnorm(n=n, sd=sd.nz), ncol=1)
  Y <- X %*% alpha + gamma * Z + rnorm(n=n, sd=sd.ny)
  
  return(list(X, Y))
}