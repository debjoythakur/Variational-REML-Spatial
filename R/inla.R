###############################################################
## INLA implementation of the Gaussian ICAR model
###############################################################

fit_inla_icar_train <- function(
    y,
    X,
    Rmat,
    train_idx,
    test_idx
) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required. Install it before running this method."
    )
  }

  n <- length(y)
  y_inla <- y
  y_inla[test_idx] <- NA

  dat_inla <- data.frame(
    y = y_inla,
    x1 = X[, 2],
    x2 = X[, 3],
    id = seq_len(n)
  )

  Cmat <- Matrix::Matrix(Rmat, sparse = TRUE)

  formula_inla <- y ~ 1 + x1 + x2 +
    f(
      id,
      model = "generic0",
      Cmatrix = Cmat,
      constr = TRUE,
      hyper = list(
        prec = list(
          prior = "loggamma",
          param = c(1, 0.01)
        )
      )
    )

  t0 <- proc.time()[3]

  fit <- INLA::inla(
    formula_inla,
    data = dat_inla,
    family = "gaussian",
    control.predictor = list(compute = TRUE),
    control.compute = list(
      config = TRUE,
      dic = FALSE,
      waic = FALSE,
      cpo = FALSE
    ),
    num.threads = 1,
    verbose = FALSE
  )

  runtime <- proc.time()[3] - t0

  prec_eps_hat <- extract_inla_hyper_mean(
    fit,
    "Precision for the Gaussian observations"
  )
  prec_u_hat <- extract_inla_hyper_mean(
    fit,
    "Precision for id"
  )

  sigma2_e_hat <- ifelse(
    is.na(prec_eps_hat),
    NA_real_,
    1 / prec_eps_hat
  )
  sigma2_u_hat <- ifelse(
    is.na(prec_u_hat),
    NA_real_,
    1 / prec_u_hat
  )

  beta_hat <- c(
    fit$summary.fixed["(Intercept)", "mean"],
    fit$summary.fixed["x1", "mean"],
    fit$summary.fixed["x2", "mean"]
  )

  uvar_all <- if ("sd" %in% colnames(fit$summary.random$id)) {
    fit$summary.random$id$sd^2
  } else {
    rep(NA_real_, n)
  }

  list(
    method = "INLA",
    sigma2_e = sigma2_e_hat,
    sigma2_u = sigma2_u_hat,
    beta = beta_hat,
    uhat_all = as.vector(fit$summary.random$id$mean),
    yhat_all = as.vector(fit$summary.fitted.values$mean),
    uvar_all = as.vector(uvar_all),
    runtime = runtime,
    converged = TRUE,
    fit_object = fit
  )
}
