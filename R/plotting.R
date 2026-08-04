###############################################################
## Simulation summaries, tables, and figures
###############################################################

summarise_simulation_results <- function(
    results_all,
    beta_true = c(1.0, 1.2, -1.0),
    sigma2_e_true = 0.7
) {
  results_all |>
    dplyr::group_by(method, m, n) |>
    dplyr::summarise(
      avg_MSPE = mean(MSPE, na.rm = TRUE),
      avg_u_mean_MSE = mean(u_mean_MSPE, na.rm = TRUE),
      avg_u_variance_MSE = mean(u_var_RMSE^2, na.rm = TRUE),
      avg_beta_MSE = mean(
        (beta0_hat - beta_true[1])^2 +
          (beta1_hat - beta_true[2])^2 +
          (beta2_hat - beta_true[3])^2,
        na.rm = TRUE
      ),
      avg_KL_q_p = mean(KL_q_p_hat, na.rm = TRUE),
      RMSE_sigma2_e = sqrt(
        mean(
          (sigma2_e_hat - sigma2_e_true)^2,
          na.rm = TRUE
        )
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(method, n)
}

make_kl_plot <- function(summary_by_n) {
  kl_plot_df <- summary_by_n |>
    dplyr::filter(method == "VRMLE") |>
    dplyr::arrange(n)

  ggplot2::ggplot(
    kl_plot_df,
    ggplot2::aes(x = n, y = avg_KL_q_p)
  ) +
    ggplot2::geom_line(
      linewidth = 0.9,
      color = "black"
    ) +
    ggplot2::geom_point(
      size = 2.4,
      color = "black"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      x = "Sample size n",
      y = expression(
        mean~KL*
          "{"*
          q(u)~"||"~g[hat(theta)](u~"|"~Y)*
          "}"
      )
    )
}

make_simulation_performance_table <- function(summary_by_n) {
  summary_by_n |>
    dplyr::transmute(
      Method = dplyr::recode(
        method,
        "Low-rank RBF VREML" = "LRCVRMLE"
      ),
      n = n,
      `Avg. MSPE` = avg_MSPE,
      `Avg. Posterior u Mean MSE` = avg_u_mean_MSE,
      `Avg. Posterior u Variance MSE` = avg_u_variance_MSE,
      `Avg. Beta MSE` = avg_beta_MSE
    ) |>
    dplyr::mutate(
      dplyr::across(
        where(is.numeric) & !all_of("n"),
        ~ round(.x, 3)
      )
    ) |>
    dplyr::arrange(
      factor(
        Method,
        levels = c(
          "INLA",
          "LRCVRMLE",
          "MLE",
          "VRMLE"
        )
      ),
      n
    )
}

make_sigma2_e_plot <- function(summary_by_n) {
  bw_linetypes <- c(
    "INLA" = "dotted",
    "MLE" = "dashed",
    "VRMLE" = "solid",
    "Low-rank RBF VREML" = "dotdash"
  )

  bw_shapes <- c(
    "INLA" = 16,
    "MLE" = 17,
    "VRMLE" = 15,
    "Low-rank RBF VREML" = 18
  )

  ggplot2::ggplot(
    summary_by_n,
    ggplot2::aes(
      x = n,
      y = RMSE_sigma2_e,
      linetype = method,
      shape = method,
      group = method
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.9,
      color = "black"
    ) +
    ggplot2::geom_point(
      size = 2.4,
      color = "black"
    ) +
    ggplot2::scale_linetype_manual(
      values = bw_linetypes,
      labels = c(
        "INLA" = "INLA",
        "MLE" = "MLE",
        "VRMLE" = "VRMLE",
        "Low-rank RBF VREML" = "LRCVRMLE"
      )
    ) +
    ggplot2::scale_shape_manual(
      values = bw_shapes,
      labels = c(
        "INLA" = "INLA",
        "MLE" = "MLE",
        "VRMLE" = "VRMLE",
        "Low-rank RBF VREML" = "LRCVRMLE"
      )
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      x = "Sample size n",
      y = expression(
        RMSE*"("*hat(sigma)[epsilon]^2*")"
      ),
      linetype = "Method",
      shape = "Method"
    )
}
