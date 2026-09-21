################################################################################
# Whole-plot design factor
################################################################################

# Internal helper for the known factor eta_{N,a} in the equal-covariance
# procedure. The calculation is kept separate so that it can be tested
# independently from the stochastic trace estimators.
compute_eta_Na <- function(TW, group_sizes, cov.equal) {
  if(!cov.equal) return(1)
  
  a <- nrow(TW)
  N <- sum(group_sizes)
  
  D_N <- diag(N / group_sizes, nrow = a, ncol = a)
  DNTW <- D_N %*% TW

  sum(diag(DNTW %*% DNTW))^3 / sum(diag(DNTW %*% DNTW %*% DNTW))^2
}
