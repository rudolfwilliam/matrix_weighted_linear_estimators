# rosenman weight

rosenman_weight <- function(alpha.hat.O, alpha.hat.I, cov.est) {
  weight <- max(c(0, 1 - (sum(diag(cov.est))/sum((alpha.hat.O$coefficients - alpha.hat.I$coefficients)^2))))
  return(weight)
}