###############################################################
## Numerical and statistical utility functions
###############################################################

safe_solve <- function(A, b = NULL, ridge = 1e-6) {
  A <- as.matrix(A)
  A <- (A + t(A)) / 2

  ridges <- c(ridge, 1e-5, 1e-4, 1e-3, 1e-2)

  for (rr in ridges) {
    A2 <- A + rr * diag(nrow(A))
    out <- tryCatch(
      {
        if (is.null(b)) solve(A2) else solve(A2, b)
      },
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }

  A2 <- A + 1e-2 * diag(nrow(A))
  if (is.null(b)) MASS::ginv(A2) else MASS::ginv(A2) %*% b
}

solve_chol <- function(A, b = NULL, ridge = 1e-6) {
  A <- as.matrix(A)
  A <- (A + t(A)) / 2

  ridges <- c(ridge, 1e-5, 1e-4, 1e-3, 1e-2)

  for (rr in ridges) {
    A2 <- A + rr * diag(nrow(A))
    out <- tryCatch(
      {
        U <- chol(A2)
        if (is.null(b)) {
          chol2inv(U)
        } else {
          backsolve(U, forwardsolve(t(U), b))
        }
      },
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }

  safe_solve(A, b, ridge = 1e-2)
}

logdet_chol <- function(A, ridge = 1e-6) {
  A <- as.matrix(A)
  A <- (A + t(A)) / 2

  ridges <- c(ridge, 1e-5, 1e-4, 1e-3, 1e-2)

  for (rr in ridges) {
    out <- tryCatch(
      {
        U <- chol(A + rr * diag(nrow(A)))
        2 * sum(log(diag(U)))
      },
      error = function(e) NULL
    )
    if (!is.null(out) && is.finite(out)) return(out)
  }

  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 1e-8)
  sum(log(ev))
}

extract_inla_hyper_mean <- function(inla_obj, pattern) {
  rn <- rownames(inla_obj$summary.hyperpar)
  idx <- grep(pattern, rn, ignore.case = TRUE)

  if (length(idx) == 0) return(NA_real_)
  inla_obj$summary.hyperpar[idx[1], "mean"]
}

kl_gaussian <- function(mu_q, Sigma_q, mu_p, Sigma_p) {
  q <- length(mu_q)
  Prec_p <- solve_chol(Sigma_p)

  term_trace <- sum(diag(Prec_p %*% Sigma_q))
  diff <- matrix(mu_p - mu_q, ncol = 1)
  term_quad <- as.numeric(t(diff) %*% Prec_p %*% diff)

  logdet_q <- logdet_chol(Sigma_q)
  logdet_p <- logdet_chol(Sigma_p)

  0.5 * (term_trace + term_quad - q + logdet_p - logdet_q)
}
