###############################################################
## Variational restricted maximum-likelihood estimation
###############################################################

fit_variational_reml_train <- function(
    y_tr,
    X_tr,
    H_tr,
    K,
    tol = 1e-8,
    max_iter = 1000
) {
  ntr <- length(y_tr)
  p <- ncol(X_tr)
  q <- ncol(H_tr)

  P_perp <- diag(ntr) -
    X_tr %*% solve_chol(t(X_tr) %*% X_tr) %*% t(X_tr)

  HPH <- t(H_tr) %*% P_perp %*% H_tr
  HPy <- t(H_tr) %*% P_perp %*% y_tr

  tau_y <- 1 / max(var(y_tr) * 0.5, 1e-3)
  tau_u <- 1 / max(var(y_tr) * 0.5, 1e-3)

  mu_theta <- rep(0, q)
  Sigma_theta <- diag(q)
  elbo_trace <- numeric(max_iter)

  converged <- FALSE
  t0 <- proc.time()[3]

  for (iter in seq_len(max_iter)) {
    Prec_q <- tau_y * HPH + tau_u * K
    Sigma_theta <- safe_solve(Prec_q)
    mu_theta <- as.vector(Sigma_theta %*% (tau_y * HPy))

    resid_mean <- y_tr - H_tr %*% mu_theta

    quad_y <- as.numeric(
      t(resid_mean) %*% P_perp %*% resid_mean
    )
    trace_y <- sum(diag(HPH %*% Sigma_theta))

    quad_u <- as.numeric(
      t(mu_theta) %*% K %*% mu_theta
    )
    trace_u <- sum(diag(K %*% Sigma_theta))

    tau_y_new <- (ntr - p) / max(quad_y + trace_y, 1e-12)
    tau_u_new <- q / max(quad_u + trace_u, 1e-12)

    logdetSigma <- logdet_chol(Sigma_theta)

    elbo_trace[iter] <-
      ((ntr - p) / 2) * log(tau_y_new) -
      (tau_y_new / 2) * (quad_y + trace_y) +
      (q / 2) * log(tau_u_new) -
      (tau_u_new / 2) * (quad_u + trace_u) +
      0.5 * logdetSigma

    if (
      iter > 1 &&
      abs(elbo_trace[iter] - elbo_trace[iter - 1]) < tol
    ) {
      tau_y <- tau_y_new
      tau_u <- tau_u_new
      elbo_trace <- elbo_trace[seq_len(iter)]
      converged <- TRUE
      break
    }

    tau_y <- tau_y_new
    tau_u <- tau_u_new

    if (iter == max_iter) {
      elbo_trace <- elbo_trace[seq_len(iter)]
    }
  }

  runtime <- proc.time()[3] - t0

  sigma2_e_hat <- 1 / tau_y
  sigma2_u_hat <- 1 / tau_u

  G_tr <- H_tr %*% solve_chol(K) %*% t(H_tr)
  V_tr <- sigma2_e_hat * diag(ntr) + sigma2_u_hat * G_tr
  Vinv_tr <- solve_chol(V_tr)

  beta_hat <- solve_chol(
    t(X_tr) %*% Vinv_tr %*% X_tr,
    t(X_tr) %*% Vinv_tr %*% y_tr
  )

  list(
    method = "VRMLE",
    sigma2_e = sigma2_e_hat,
    sigma2_u = sigma2_u_hat,
    beta = as.vector(beta_hat),
    mu_theta = mu_theta,
    Sigma_theta = Sigma_theta,
    runtime = runtime,
    converged = converged,
    n_iter = length(elbo_trace),
    elbo_trace = elbo_trace
  )
}
