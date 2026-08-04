###############################################################
## Main simulation driver
##
## Run this script from the repository root:
## source("simulations/run_simulation.R")
###############################################################

rm(list = ls())

required_packages <- c(
  "MASS",
  "Matrix",
  "ggplot2",
  "dplyr",
  "RANN"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the simulation: ",
    paste(missing_packages, collapse = ", ")
  )
}

if (!requireNamespace("INLA", quietly = TRUE)) {
  stop(
    "Package 'INLA' is required. Install it from the official INLA repository."
  )
}

source("simulations/simulation_parameters.R")

source("R/utilities.R")
source("R/graph_generation.R")
source("R/data_generation.R")
source("R/exact_mle.R")
source("R/vrmle.R")
source("R/lowrank_vrmle.R")
source("R/inla.R")
source("R/posterior.R")
source("R/prediction.R")
source("R/plotting.R")

set.seed(simulation_config$seed)

dir.create(
  simulation_config$output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(simulation_config$output_dir, "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(simulation_config$output_dir, "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(simulation_config$output_dir, "output"),
  recursive = TRUE,
  showWarnings = FALSE
)

results_list <- list()
result_counter <- 0L

for (m in simulation_config$m_values) {
  cat("\n============================================\n")
  cat("Running simulations for m =", m, "=> n =", m^2, "\n")
  cat("============================================\n")

  scn <- build_scenario(
    m = m,
    beta_true = simulation_config$beta_true,
    sigma2_e_true = simulation_config$sigma2_e_true,
    sigma2_u_true = simulation_config$sigma2_u_true
  )

  for (s in seq_len(simulation_config$nsim)) {
    cat(
      "Simulation",
      s,
      "of",
      simulation_config$nsim,
      "for n =",
      scn$n,
      "\n"
    )

    tmp <- tryCatch(
      one_run_prediction(
        scn = scn,
        train_frac = simulation_config$train_frac,
        k_basis = simulation_config$k_basis,
        seed = simulation_config$seed + 1000L * m + s
      ),
      error = function(e) {
        message(
          "Simulation failed for m=",
          m,
          ", sim=",
          s,
          ": ",
          e$message
        )
        NULL
      }
    )

    if (!is.null(tmp)) {
      tmp$sim <- s
      result_counter <- result_counter + 1L
      results_list[[result_counter]] <- tmp
    }
  }
}

if (length(results_list) == 0) {
  stop("All simulation runs failed; no results were generated.")
}

results_all <- dplyr::bind_rows(results_list)

summary_by_n <- summarise_simulation_results(
  results_all = results_all,
  beta_true = simulation_config$beta_true,
  sigma2_e_true = simulation_config$sigma2_e_true
)

simulation_performance <- make_simulation_performance_table(
  summary_by_n
)

figure_1_KL <- make_kl_plot(summary_by_n)
figure_2_sigma2_e <- make_sigma2_e_plot(summary_by_n)

print(simulation_performance)
print(figure_1_KL)
print(figure_2_sigma2_e)

write.csv(
  results_all,
  file.path(
    simulation_config$output_dir,
    "output",
    "simulation_replication_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_by_n,
  file.path(
    simulation_config$output_dir,
    "tables",
    "simulation_summary_by_n.csv"
  ),
  row.names = FALSE
)

write.csv(
  simulation_performance,
  file.path(
    simulation_config$output_dir,
    "tables",
    "simulation_performance_table.csv"
  ),
  row.names = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    simulation_config$output_dir,
    "figures",
    "figure_1_KL.pdf"
  ),
  plot = figure_1_KL,
  width = 7,
  height = 5
)

ggplot2::ggsave(
  filename = file.path(
    simulation_config$output_dir,
    "figures",
    "figure_2_sigma2_e.pdf"
  ),
  plot = figure_2_sigma2_e,
  width = 7,
  height = 5
)

saveRDS(
  list(
    results_all = results_all,
    summary_by_n = summary_by_n,
    simulation_performance = simulation_performance,
    configuration = simulation_config
  ),
  file.path(
    simulation_config$output_dir,
    "output",
    "simulation_results.rds"
  )
)

cat(
  "\nSimulation completed. Outputs were saved under:",
  simulation_config$output_dir,
  "\n"
)
