#' Niche analysis using ECDF and chi-squared
#'
#' Simulate niche suitability from Mahalanobis distance using both
#' chi-squared and empirical CDF transformations, for a given number
#' of predictor variables.
#'
#' @param n Integer; number of predictor variables (dimensions).
#' @param n_population Integer; size of simulated environmental population.
#' @param sample_sizes Integer vector of sample sizes to evaluate.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#' \itemize{
#'   \item corplot: ggplot object with correlation vs sample size.
#'   \item sample_data: matrix of simulated sample points.
#'   \item sample_niche: numeric vector of “true” niche suitability.
#'   \item chisq_suits: numeric vector, 1 - pchisq(Mahalanobis).
#'   \item ecdf_suits: numeric vector, 1 - ECDF(Mahalanobis).
#'   \item mahal_dists: numeric vector of Mahalanobis distances.
#' }
#'
#' @author Luíz Fernando Esser
#'
#' @examples
#' # Create ECDF-niche based on personalized options:
#' n <- ecdf_theoretical_niche(n = 3,
#'                             n_population = 20000,
#'                             sample_sizes = seq(50, 1000, 50),
#'                             seed = 123)
#'
#' @importFrom MASS mvrnorm
#' @importFrom ggplot2 ggplot aes geom_line theme_bw labs scale_colour_manual theme element_blank
#'             ylim
#' @importFrom stats pchisq ecdf cor mahalanobis cov
#' @import checkCLI
#'
#' @global .data
#'
#' @export
ecdf_theoretical_niche <- function(
    n,
    n_population = 10000L,
    sample_sizes = seq(20L, 500L, 20L),
    seed = NULL) {

  # Assertions
  assert_numeric_cli(n, len = 1, lower = 1)
  assert_numeric_cli(n_population, len = 1, lower = 1)
  assert_vector_cli(sample_sizes, min.len = 1, unique = TRUE)
  assert_numeric_cli(sample_sizes, lower = 1L, null.ok = TRUE)
  assert_numeric_cli(seed, lower = 0, null.ok = TRUE)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # population
  population_grid <- lapply(seq_len(n), function(x) stats::rnorm(n_population))
  population_grid <- as.data.frame(do.call(cbind, population_grid))
  population_means <- colMeans(population_grid)
  population_covariance <- stats::cov(population_grid)

  ecdf_correlations  <- numeric(0L)
  chisq_correlations <- numeric(0L)

  # store last sample for returning
  last_sample_data         <- NULL
  last_niche_suitability   <- NULL
  last_chisq_suitabilities <- NULL
  last_ecdf_suitabilities  <- NULL
  last_mahalanobis         <- NULL

  for (current_sample_size in sample_sizes) {
    sample_data <- MASS::mvrnorm(
      n   = current_sample_size,
      mu  = population_means,
      Sigma = population_covariance
    )

    mahalanobis_distances <- stats::mahalanobis(
      x      = sample_data,
      center = population_means,
      cov    = population_covariance
    )

    niche_suitability <- 1 - stats::pchisq(mahalanobis_distances, df = n)

    chisq_suitabilities <- 1 - stats::pchisq(mahalanobis_distances, df = n)

    ecdf_function <- stats::ecdf(mahalanobis_distances)
    ecdf_suitabilities <- 1 - ecdf_function(mahalanobis_distances)

    correlation_matrix <- stats::cor(cbind(
      niche_suitability,
      chisq_suitabilities,
      ecdf_suitabilities
    ))

    ecdf_correlations  <- c(ecdf_correlations,  correlation_matrix[1, 3])
    chisq_correlations <- c(chisq_correlations, correlation_matrix[1, 2])

    last_sample_data         <- sample_data
    last_niche_suitability   <- niche_suitability
    last_chisq_suitabilities <- chisq_suitabilities
    last_ecdf_suitabilities  <- ecdf_suitabilities
    last_mahalanobis         <- mahalanobis_distances
  }

  plot_data <- data.frame(
    Sample_Size = c(sample_sizes, sample_sizes),
    Correlation = c(chisq_correlations, ecdf_correlations),
    Method      = c(
      rep("1 - Chi-squared", length(sample_sizes)),
      rep("1 - ECDF",       length(sample_sizes))
    )
  )

  correlation_plot <- ggplot2::ggplot(data = plot_data) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
      ggplot2::aes(x = .data$Sample_Size, y = .data$Correlation, colour = .data$Method),
      linewidth = 1.0,
      stat = "identity"
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Number of Records",
      y = "Correlation",
      title = paste(n, "Predictor Variables")
    ) +
    ggplot2::scale_colour_manual(values = c("cyan4", "orange"))

  list(
    corplot     = correlation_plot,
    sample_data = last_sample_data,
    sample_niche = last_niche_suitability,
    chisq_suits  = last_chisq_suitabilities,
    ecdf_suits   = last_ecdf_suitabilities,
    mahal_dists  = last_mahalanobis
  )
}
