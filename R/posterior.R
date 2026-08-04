###############################################################
## Exact Gaussian posterior calculations
###############################################################

exact_posterior_theta <- function(
    y_tr,
    X_tr,
    H_tr,
    K,
    sigma2_e,
    sigma2_u
) {
  ntr <- length(y_tr)

  P_perp <- diag(ntr) -
    X_tr %*% solve_chol(t(X_tr) %*% X_tr) %*% t(X_tr)

  tau_y <- 1 / sigma2_e
  tau_u <- 1 / sigma2_u

  A <- tau_y * t(H_tr) %*% P_perp %*% H_tr +
    tau_u * K
  b <- tau_y * t(H_tr) %*% P_perp %*% y_tr

  Sigma <- solve_chol(A)
  mu <- as.vector(Sigma %*% b)

  list(mu = mu, Sigma = Sigma)
}
