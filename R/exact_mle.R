###############################################################
## Exact Gaussian maximum-likelihood estimation
###############################################################

fit_exact_gaussian_train <- function(y_tr, X_tr, G_tr) {
  ntr <- length(y_tr)

  obj_fun <- function(log_par) {
    sigma2_e <- exp(log_par[1])
    sigma2_u <- exp(log_par[2])

    V <- sigma2_e * diag(ntr) + sigma2_u * G_tr
    Vinv <- solve_chol(V)

    XtVinvX <- t(X_tr) %*% Vinv %*% X_tr
    beta_hat <- solve_chol(
      XtVinvX,
      t(X_tr) %*% Vinv %*% y_tr
    )

    resid <- y_tr - X_tr %*% beta_hat

    ll <- -0.5 * (
      ntr * log(2 * pi) +
        logdet_chol(V) +
        as.numeric(t(resid) %*% Vinv %*% resid)
    )

    -ll
  }

  vy <- var(y_tr)
  init <- log(c(
    max(vy * 0.4, 1e-3),
    max(vy * 0.6, 1e-3)
  ))

  t0 <- proc.time()[3]

  opt <- optim(
    par = init,
    fn = obj_fun,
    method = "L-BFGS-B",
    lower = log(c(1e-6, 1e-6)),
    upper = log(c(1e3, 1e3))
  )

  runtime <- proc.time()[3] - t0

  sigma2_e_hat <- exp(opt$par[1])
  sigma2_u_hat <- exp(opt$par[2])

  V <- sigma2_e_hat * diag(ntr) + sigma2_u_hat * G_tr
  Vinv <- solve_chol(V)

  beta_hat <- solve_chol(
    t(X_tr) %*% Vinv %*% X_tr,
    t(X_tr) %*% Vinv %*% y_tr
  )

  list(
    method = "MLE",
    sigma2_e = sigma2_e_hat,
    sigma2_u = sigma2_u_hat,
    beta = as.vector(beta_hat),
    runtime = runtime,
    converged = opt$convergence == 0
  )
}
