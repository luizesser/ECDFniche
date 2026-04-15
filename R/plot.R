#' Create distance–suitability plot
#'
#' @param analysis_results List returned by \code{ecdf_theoretical_niche()}.
#' @return A ggplot object.
#'
#' @examples
#' # Create ECDF-niche based on personalized options:
#' res <- ecdf_theoretical_niche(n = 3,
#'                               n_population = 20000,
#'                               sample_sizes = seq(50, 1000, 50),
#'                               seed = 123)
#'
#' # Plot analysis results
#' create_distance_suitability_plot(res)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_colour_manual theme_classic labs theme element_text
#' @importFrom stats pchisq ecdf
#' @import checkCLI
#'
#' @global .data
#'
#' @export
create_distance_suitability_plot <- function(analysis_results) {
  # Assertions
  assert_list_cli(analysis_results, len = 6)
  assert_names_cli(names(analysis_results), must.include = c("corplot", "sample_data",
                                                      "sample_niche", "chisq_suits",
                                                      "ecdf_suits", "mahal_dists"))
  assert_class_cli(analysis_results$corplot, classes = "ggplot")
  assert_matrix_cli(analysis_results$sample_data)
  assert_numeric_cli(analysis_results$sample_niche, lower = 0, upper = 1)
  assert_numeric_cli(analysis_results$chisq_suits, lower = 0, upper = 1)
  assert_numeric_cli(analysis_results$ecdf_suits, lower = 0, upper = 1)
  assert_numeric_cli(analysis_results$mahal_dists, lower = 0)

  plot_data <- data.frame(
    Mahalanobis_Distance   = analysis_results$mahal_dists,
    Niche_Suitability      = analysis_results$sample_niche,
    ChiSquared_suitability = analysis_results$chisq_suits,
    ECDF_suitability       = analysis_results$ecdf_suits
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Mahalanobis_Distance)) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$Niche_Suitability, color = "Niche Records"),
      shape = 21, size = 3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$ChiSquared_suitability, color = "1 - Chi-squared"),
      linewidth = 1
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$ECDF_suitability, color = "1 - ECDF"),
      linewidth = 1
    ) +
    ggplot2::scale_colour_manual(
      values = c("Niche Records" = "black",
                 "1 - Chi-squared" = "cyan4",
                 "1 - ECDF" = "orange")
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Mahalanobis Distance",
      y = "Environmental Suitability",
      title = paste0(ncol(analysis_results$sample_data), " Predictor Variables")
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme(
      axis.title  = ggplot2::element_text(size = 12),
      axis.text   = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_blank(),
      legend.position = c(0.8, 0.8)
    )
}

#' Run full ECDF–Mahalanobis analysis
#'
#' Convenience function that reproduces the three figures from the original
#' manuscript for 1–5 dimensions.
#'
#' @param dims Integer vector of dimensions (default 1:5).
#' @param seed Optional seed for reproducibility.
#'
#' @return A list containing:
#' \itemize{
#'   \item analyses: list of ecdf_theoretical_niche() outputs.
#'   \item figure1, figure2, figure3: grobs with arranged plots.
#' }
#'
#' @examples
#' # Recreate original manuscript output:
#' set.seed(3)
#' full_res <- run_ecdf_mahal_analysis(dims = 1:5)
#'
#' @importFrom ggplot2 ggplot aes geom_point labs scale_colour_viridis_c theme_classic theme element_blank
#' @importFrom lemon grid_arrange_shared_legend
#' @import checkCLI
#'
#' @global .data
#'
#' @export
run_ecdf_mahal_analysis <- function(dims = 1:5, seed = 3L) {
  # Assertions
  assert_vector_cli(dims, min.len = 1, unique = TRUE)
  assert_numeric_cli(dims, lower = 1)
  assert_numeric_cli(seed, lower = 0, null.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  analyses <- lapply(dims, function(d) ecdf_theoretical_niche(d, seed = seed))

  vars1 <- analyses[[1L]]
  vars2 <- analyses[[2L]]
  vars3 <- analyses[[3L]]
  vars4 <- analyses[[4L]]
  vars5 <- analyses[[5L]]

  # Figure 1 (2D spatial)
  spatial_data <- data.frame(
    Variable_1             = vars2$sample_data[, 1],
    Variable_2             = vars2$sample_data[, 2],
    Niche_suitability      = vars2$sample_niche,
    ChiSquared_suitability = vars2$chisq_suits,
    ECDF_suitability       = vars2$ecdf_suits
  )

  common_theme <- ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none")

  plot_niche <- ggplot2::ggplot(spatial_data,
                                ggplot2::aes(x = .data$Variable_1, y = .data$Variable_2)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Niche_suitability), size = 3) +
    ggplot2::labs(
      x = "",
      y = "Variable 2",
      title = "Simulated Niche",
      color = "Environmental Suitability"
    ) +
    ggplot2::scale_colour_viridis_c() +
    common_theme

  plot_chisq <- ggplot2::ggplot(spatial_data,
                                ggplot2::aes(x = .data$Variable_1, y = .data$Variable_2)) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$ChiSquared_suitability), size = 3
    ) +
    ggplot2::labs(
      x = "Variable 1",
      y = "",
      title = "1 - Chi-squared Suitability"
    ) +
    ggplot2::scale_colour_viridis_c() +
    common_theme

  plot_ecdf <- ggplot2::ggplot(spatial_data,
                               ggplot2::aes(x = .data$Variable_1, y = .data$Variable_2)) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$ECDF_suitability), size = 3
    ) +
    ggplot2::labs(
      x = "",
      y = "",
      title = "1 - ECDF Suitability"
    ) +
    ggplot2::scale_colour_viridis_c() +
    common_theme

  figure1 <- lemon::grid_arrange_shared_legend(
    plot_niche, plot_chisq, plot_ecdf,
    ncol = 3, nrow = 1, position = "bottom"
  )

  # Figure 2: correlation plots
  figure2 <- lemon::grid_arrange_shared_legend(
    vars1$corplot + ggplot2::labs(x = "", y = "Correlation") +
      ggplot2::coord_cartesian(ylim = c(0.95, 1)),
    vars2$corplot + ggplot2::labs(x = "", y = "") +
      ggplot2::coord_cartesian(ylim = c(0.95, 1)),
    vars3$corplot + ggplot2::labs(x = "Number of Records", y = "") +
      ggplot2::coord_cartesian(ylim = c(0.95, 1)),
    vars4$corplot + ggplot2::labs(x = "", y = "") +
      ggplot2::coord_cartesian(ylim = c(0.95, 1)),
    vars5$corplot + ggplot2::labs(x = "", y = "") +
      ggplot2::coord_cartesian(ylim = c(0.95, 1)),
    ncol = 5, nrow = 1, position = "bottom"
  )

  # Figure 3: distance–suitability relationships
  distance_plot_1 <- create_distance_suitability_plot(vars1)
  distance_plot_2 <- create_distance_suitability_plot(vars2)
  distance_plot_3 <- create_distance_suitability_plot(vars3)
  distance_plot_4 <- create_distance_suitability_plot(vars4)
  distance_plot_5 <- create_distance_suitability_plot(vars5)

  figure3 <- lemon::grid_arrange_shared_legend(
    distance_plot_1 + ggplot2::labs(x = "", y = "Environmental Suitability"),
    distance_plot_2 + ggplot2::labs(x = "", y = ""),
    distance_plot_3 + ggplot2::labs(x = "Mahalanobis Distance", y = ""),
    distance_plot_4 + ggplot2::labs(x = "", y = ""),
    distance_plot_5 + ggplot2::labs(x = "", y = ""),
    ncol = 5, nrow = 1, position = "bottom"
  )

  list(
    analyses = analyses,
    figure1  = figure1,
    figure2  = figure2,
    figure3  = figure3
  )
}
