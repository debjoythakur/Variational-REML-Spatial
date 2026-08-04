###############################################################
## Gaussian ICAR data generation
###############################################################

simulate_icar_data <- function(X, H, Kinv, beta, sigma2_e, sigma2_u) {
  n <- nrow(X)
  q <- ncol(H)

  theta <- as.vector(
    MASS::mvrnorm(
      n = 1,
      mu = rep(0, q),
      Sigma = sigma2_u * Kinv
    )
  )

  u <- as.vector(H %*% theta)
  eps <- rnorm(n, mean = 0, sd = sqrt(sigma2_e))
  y <- as.vector(X %*% beta + u + eps)

  list(y = y, u = u, theta = theta)
}
