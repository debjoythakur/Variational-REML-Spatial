###############################################################
## User-adjustable simulation parameters
###############################################################

simulation_config <- list(
  seed = 123,

  ## Grid dimensions; sample size is n = m^2
  m_values = c(5, 8, 10, 20, 25, 30),

  ## Number of Monte Carlo replications
  nsim = 10,

  ## Proportion of observations used for model fitting
  train_frac = 0.70,

  ## Number of low-rank RBF basis functions
  k_basis = 20,

  ## Data-generating parameters
  beta_true = c(1.0, 1.2, -1.0),
  sigma2_e_true = 0.7,
  sigma2_u_true = 1.3,

  ## Output directory
  output_dir = "results"
)
