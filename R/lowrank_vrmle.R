###############################################################
## Low-rank radial-basis-function VREML
###############################################################

make_rbf_basis <- function(
    coords,
    k_basis = 20,
    lengthscale = NULL,
    seed = 123
) {
  set.seed(seed)

  coords <- as.matrix(coords)
  n <- nrow(coords)
  k_basis <- min(k_basis, max(2, floor(n / 3)))

  km <- tryCatch(
    kmeans(
      coords,
      centers = k_basis,
      nstart = 25,
      iter.max = 100
    ),
    warning = function(w) {
      kmeans(
        coords,
        centers = k_basis,
        nstart = 25,
        iter.max = 100
      )
    }
  )

  centers <- as.matrix(km$centers)

  dist_xc <- as.matrix(dist(rbind(coords, centers)))
  D_xc <- dist_xc[
    seq_len(n),
    n + seq_len(k_basis),
    drop = FALSE
  ]

  if (is.null(lengthscale)) {
    cc <- as.matrix(dist(centers))
    lengthscale <- median(cc[upper.tri(cc)], na.rm = TRUE)

    if (!is.finite(lengthscale) || lengthscale <= 0) {
      lengthscale <- median(D_xc, na.rm = TRUE)
    }
  }

  Phi_raw <- exp(-(D_xc^2) / (2 * lengthscale^2))
  Phi_mean <- colMeans(Phi_raw)
  Phi <- sweep(Phi_raw, 2, Phi_mean, "-")

  keep <- apply(Phi, 2, sd) > 1e-8
  Phi <- Phi[, keep, drop = FALSE]
  Phi_mean <- Phi_mean[keep]
  centers <- centers[keep, , drop = FALSE]

  qrP <- qr(Phi)
  Phi_orth <- qr.Q(qrP)[, seq_len(qrP$rank), drop = FALSE]

  list(
    Phi = Phi_orth,
    centers = centers,
    lengthscale = lengthscale,
    Phi_mean = Phi_mean,
    original_Phi = Phi
  )
}

fit_lowrank_rbf_vreml_train <- function(
    y_tr,
    X_tr,
    coords_tr,
    R_tr,
    k_basis = 20,
    lengthscale = NULL,
    tol = 1e-8,
    max_iter = 1000,
    seed = 123
) {
  ntr <- length(y_tr)
  p <- ncol(X_tr)

  k_basis_use <- min(
    k_basis,
    floor(ntr / 3),
    ntr - p - 2
  )
  k_basis_use <- max(2, k_basis_use)

  rbf <- make_rbf_basis(
    coords = coords_tr,
    k_basis = k_basis_use,
    lengthscale = lengthscale,
    seed = seed
  )

  Phi <- rbf$Phi
  k <- ncol(Phi)

  P_perp <- diag(ntr) -
    X_tr %*% solve_chol(t(X_tr) %*% X_tr) %*% t(X_tr)

  PhiP <- t(Phi) %*% P_perp %*% Phi
  PhiR <- t(Phi) %*% R_tr %*% Phi
  PhiPy <- t(Phi) %*% P_perp %*% y_tr

  tau_y <- 1 / max(var(y_tr) * 0.5, 1e-3)
  tau_u <- 1 / max(var(y_tr) * 0.5, 1e-3)

  mu_a <- rep(0, k)
  Sigma_a <- diag(k)
  elbo_trace <- numeric(max_iter)

  converged <- FALSE
  t0 <- proc.time()[3]

  for (iter in seq_len(max_iter)) {
    Prec_a <- tau_y * PhiP + tau_u * PhiR
    Sigma_a <- safe_solve(Prec_a)
    mu_a <- as.vector(Sigma_a %*% (tau_y * PhiPy))

    uhat <- as.vector(Phi %*% mu_a)

    quad_y <- as.numeric(
      t(y_tr - uhat) %*% P_perp %*% (y_tr - uhat)
    )
    trace_y <- sum(diag(PhiP %*% Sigma_a))

    quad_u <- as.numeric(
      t(mu_a) %*% PhiR %*% mu_a
    )
    trace_u <- sum(diag(PhiR %*% Sigma_a))

    tau_y_new <- (ntr - p) / max(quad_y + trace_y, 1e-12)
    tau_u_new <- k / max(quad_u + trace_u, 1e-12)

    logdetSigma <- logdet_chol(Sigma_a)

    elbo_trace[iter] <-
      ((ntr - p) / 2) * log(tau_y_new) -
      (tau_y_new / 2) * (quad_y + trace_y) +
      (k / 2) * log(tau_u_new) -
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

  uhat_tr <- as.vector(Phi %*% mu_a)
  beta_hat <- solve_chol(
    t(X_tr) %*% X_tr,
    t(X_tr) %*% (y_tr - uhat_tr)
  )

  list(
    method = "Low-rank RBF VREML",
    sigma2_e = sigma2_e_hat,
    sigma2_u = sigma2_u_hat,
    beta = as.vector(beta_hat),
    mu_a = mu_a,
    Sigma_a = Sigma_a,
    Phi = Phi,
    k_basis_used = k,
    runtime = runtime,
    converged = converged,
    n_iter = length(elbo_trace),
    elbo_trace = elbo_trace
  )
}
