###############################################################
## Prediction functions and one simulation replication
###############################################################

predict_exact_gaussian <- function(
    fit,
    y_tr,
    X_tr,
    X_te,
    G_tr,
    G_te_tr
) {
  sigma2_e <- fit$sigma2_e
  sigma2_u <- fit$sigma2_u
  beta_hat <- fit$beta

  V_tr <- sigma2_e * diag(length(y_tr)) + sigma2_u * G_tr
  Vinv_tr <- solve_chol(V_tr)
  resid_tr <- y_tr - X_tr %*% beta_hat

  uhat_te <- as.vector(
    sigma2_u * G_te_tr %*% Vinv_tr %*% resid_tr
  )
  yhat_te <- as.vector(X_te %*% beta_hat + uhat_te)

  G_te <- attr(G_te_tr, "G_te")

  Sigma_cond <- sigma2_u * G_te -
    sigma2_u^2 *
      G_te_tr %*% Vinv_tr %*% t(G_te_tr)

  list(
    yhat = yhat_te,
    uhat = uhat_te,
    uvar = diag(Sigma_cond)
  )
}

predict_variational_reml <- function(fit, X_te, H_te) {
  uhat_te <- as.vector(H_te %*% fit$mu_theta)
  yhat_te <- as.vector(X_te %*% fit$beta + uhat_te)

  Sigma_u_te <- H_te %*% fit$Sigma_theta %*% t(H_te)

  list(
    yhat = yhat_te,
    uhat = uhat_te,
    uvar = diag(Sigma_u_te)
  )
}

predict_lowrank_rbf_vreml <- function(
    fit,
    X_te,
    coords_te,
    coords_tr
) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop(
      "Package 'RANN' is required for low-rank nearest-neighbor prediction."
    )
  }

  nn <- RANN::nn2(
    data = coords_tr,
    query = coords_te,
    k = 1
  )$nn.idx[, 1]

  uhat_tr <- as.vector(fit$Phi %*% fit$mu_a)
  uhat_te <- uhat_tr[nn]
  yhat_te <- as.vector(X_te %*% fit$beta + uhat_te)

  Sigma_u_tr <- fit$Phi %*% fit$Sigma_a %*% t(fit$Phi)
  uvar_te <- diag(Sigma_u_tr)[nn]

  list(
    yhat = yhat_te,
    uhat = uhat_te,
    uvar = uvar_te
  )
}

predict_inla_icar <- function(fit, test_idx) {
  list(
    yhat = fit$yhat_all[test_idx],
    uhat = fit$uhat_all[test_idx],
    uvar = fit$uvar_all[test_idx]
  )
}

one_run_prediction <- function(
    scn,
    train_frac = 0.7,
    k_basis = 20,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  dat <- simulate_icar_data(
    X = scn$X,
    H = scn$H,
    Kinv = scn$Kinv,
    beta = scn$beta_true,
    sigma2_e = scn$sigma2_e_true,
    sigma2_u = scn$sigma2_u_true
  )

  y <- dat$y
  u <- dat$u
  n <- scn$n

  train_idx <- sort(
    sample(
      seq_len(n),
      floor(train_frac * n),
      replace = FALSE
    )
  )
  test_idx <- setdiff(seq_len(n), train_idx)

  y_tr <- y[train_idx]
  y_te <- y[test_idx]
  u_te <- u[test_idx]

  X_tr <- scn$X[train_idx, , drop = FALSE]
  X_te <- scn$X[test_idx, , drop = FALSE]

  H_tr <- scn$H[train_idx, , drop = FALSE]
  H_te <- scn$H[test_idx, , drop = FALSE]

  coords_tr <- scn$coords[train_idx, , drop = FALSE]
  coords_te <- scn$coords[test_idx, , drop = FALSE]
  R_tr <- scn$Rmat[train_idx, train_idx, drop = FALSE]

  G_tr <- scn$Gmat[train_idx, train_idx, drop = FALSE]
  G_te <- scn$Gmat[test_idx, test_idx, drop = FALSE]
  G_te_tr <- scn$Gmat[test_idx, train_idx, drop = FALSE]
  attr(G_te_tr, "G_te") <- G_te

  fit_vrmle <- fit_variational_reml_train(
    y_tr,
    X_tr,
    H_tr,
    scn$Kmat
  )

  fit_lr <- fit_lowrank_rbf_vreml_train(
    y_tr = y_tr,
    X_tr = X_tr,
    coords_tr = coords_tr,
    R_tr = R_tr,
    k_basis = min(
      k_basis,
      floor(length(y_tr) / 3),
      length(y_tr) - 5
    ),
    seed = 123
  )

  fit_mle <- fit_exact_gaussian_train(
    y_tr,
    X_tr,
    G_tr
  )

  fit_inla <- fit_inla_icar_train(
    y = y,
    X = scn$X,
    Rmat = scn$Rmat,
    train_idx = train_idx,
    test_idx = test_idx
  )

  pred_vrmle <- predict_variational_reml(
    fit_vrmle,
    X_te,
    H_te
  )
  pred_lr <- predict_lowrank_rbf_vreml(
    fit_lr,
    X_te,
    coords_te,
    coords_tr
  )
  pred_mle <- predict_exact_gaussian(
    fit_mle,
    y_tr,
    X_tr,
    X_te,
    G_tr,
    G_te_tr
  )
  pred_inla <- predict_inla_icar(
    fit_inla,
    test_idx
  )

  post_true <- exact_posterior_theta(
    y_tr = y_tr,
    X_tr = X_tr,
    H_tr = H_tr,
    K = scn$Kmat,
    sigma2_e = scn$sigma2_e_true,
    sigma2_u = scn$sigma2_u_true
  )

  Sigma_u_true_te <- H_te %*% post_true$Sigma %*% t(H_te)
  uvar_true_te <- diag(Sigma_u_true_te)

  post_hat_vrmle <- exact_posterior_theta(
    y_tr = y_tr,
    X_tr = X_tr,
    H_tr = H_tr,
    K = scn$Kmat,
    sigma2_e = fit_vrmle$sigma2_e,
    sigma2_u = fit_vrmle$sigma2_u
  )

  KL_q_p_hat <- kl_gaussian(
    mu_q = fit_vrmle$mu_theta,
    Sigma_q = fit_vrmle$Sigma_theta,
    mu_p = post_hat_vrmle$mu,
    Sigma_p = post_hat_vrmle$Sigma
  )

  KL_q_p_true <- kl_gaussian(
    mu_q = fit_vrmle$mu_theta,
    Sigma_q = fit_vrmle$Sigma_theta,
    mu_p = post_true$mu,
    Sigma_p = post_true$Sigma
  )

  make_row <- function(
      method_name,
      fit,
      pred,
      KL_hat = NA_real_,
      KL_true = NA_real_
  ) {
    data.frame(
      method = method_name,
      m = scn$m,
      n = scn$n,
      train_n = length(y_tr),
      test_n = length(y_te),
      sigma2_e_hat = fit$sigma2_e,
      sigma2_u_hat = fit$sigma2_u,
      beta0_hat = fit$beta[1],
      beta1_hat = fit$beta[2],
      beta2_hat = fit$beta[3],
      bias_beta0 = fit$beta[1] - scn$beta_true[1],
      bias_beta1 = fit$beta[2] - scn$beta_true[2],
      bias_beta2 = fit$beta[3] - scn$beta_true[3],
      MSPE = mean((y_te - pred$yhat)^2, na.rm = TRUE),
      MAE = mean(abs(y_te - pred$yhat), na.rm = TRUE),
      u_mean_MSPE = mean((u_te - pred$uhat)^2, na.rm = TRUE),
      u_var_RMSE = sqrt(
        mean((pred$uvar - uvar_true_te)^2, na.rm = TRUE)
      ),
      KL_q_p_hat = KL_hat,
      KL_q_p_true = KL_true,
      runtime = fit$runtime,
      converged = fit$converged,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    make_row(
      "VRMLE",
      fit_vrmle,
      pred_vrmle,
      KL_q_p_hat,
      KL_q_p_true
    ),
    make_row(
      "Low-rank RBF VREML",
      fit_lr,
      pred_lr
    ),
    make_row(
      "MLE",
      fit_mle,
      pred_mle
    ),
    make_row(
      "INLA",
      fit_inla,
      pred_inla
    )
  )
}
