###############################################################
## Graph construction and simulation-scenario preparation
###############################################################

make_grid_graph <- function(m) {
  n <- m * m
  W <- matrix(0, n, n)

  idx <- function(i, j) (i - 1) * m + j

  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      a <- idx(i, j)

      if (i > 1) W[a, idx(i - 1, j)] <- 1
      if (i < m) W[a, idx(i + 1, j)] <- 1
      if (j > 1) W[a, idx(i, j - 1)] <- 1
      if (j < m) W[a, idx(i, j + 1)] <- 1
    }
  }

  W <- (W + t(W)) / 2
  D <- diag(rowSums(W))
  R <- D - W

  list(W = W, D = D, R = R)
}

make_H_basis <- function(n) {
  Cn <- diag(n) - matrix(1, n, n) / n
  ee <- eigen(Cn, symmetric = TRUE)
  ee$vectors[, ee$values > 1e-8, drop = FALSE]
}

build_scenario <- function(
    m,
    beta_true = c(1.0, 1.2, -1.0),
    sigma2_e_true = 0.7,
    sigma2_u_true = 1.3
) {
  graph <- make_grid_graph(m)
  Rmat <- graph$R
  n <- nrow(Rmat)

  H <- make_H_basis(n)
  Kmat <- t(H) %*% Rmat %*% H
  Kinv <- solve_chol(Kmat)
  Gmat <- H %*% Kinv %*% t(H)

  coords <- expand.grid(row = seq_len(m), col = seq_len(m))
  x1 <- scale(coords$row)[, 1]
  x2 <- scale(coords$col)[, 1]
  X <- cbind(1, x1, x2)

  list(
    m = m,
    n = n,
    coords = as.matrix(coords),
    Rmat = Rmat,
    Wmat = graph$W,
    H = H,
    Kmat = Kmat,
    Kinv = Kinv,
    Gmat = Gmat,
    X = X,
    p = ncol(X),
    beta_true = beta_true,
    sigma2_e_true = sigma2_e_true,
    sigma2_u_true = sigma2_u_true
  )
}
