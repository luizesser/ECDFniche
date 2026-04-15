#' Simulations and analyses of Mahalanobis distance-based habitat suitability
#'
#' Performs replicated simulations of multivariate normal data to evaluate
#' the agreement between suitability derived from chi-squared distribution
#' and empirical cumulative distribution function (ECDF).
#'
#' @param p_vals Integer vector; number of predictor variables (dimensions).
#' @param n_vals Integer vector; number of records (sample sizes).
#' @param n_reps Integer; number of replicates per combination.
#' @param seed Optional integer for reproducibility.
#'
#' @description
#' The objective is to compare the performance of habitat suitability calculated
#' based on chi-squared cumulative distribution function and Empirical
#' Cumulative Distribution Function (ECDF)
#'
#' @return A list with:
#' \itemize{
#'   \item cor_plot: ggplot of correlation vs sample size.
#'   \item suit_plot: ggplot of suitability vs Mahalanobis distance.
#'   \item cond_plot: ggplot of correlation vs condition number.
#'   \item cor_df: raw correlation data.
#'   \item obs_df: observation-level data.
#'   \item cov_df: covariance diagnostics.
#' }
#'
#' @author Matheus T. Baumgartner
#'
#' @examples
#' # Create ECDF-niche based on personalized options:
#' n <- ecdf_compare_niche(p_vals = 1:5,
#'                         n_vals = seq(20L, 500L, 20L),
#'                         n_reps = 30L,
#'                         seed = 1991)
#'
#' @importFrom stats rWishart mahalanobis cov pchisq ecdf cor cov2cor sd
#' @importFrom MASS mvrnorm
#' @importFrom dplyr group_by summarise mutate ungroup sample_frac
#' @importFrom ggplot2 ggplot geom_point geom_pointrange stat_function
#'   facet_wrap labs theme_bw theme element_blank element_text scale_color_manual
#'   label_both
#' @import checkCLI
#'
#' @global .data
#'
#' @export
ecdf_compare_niche <- function(
    p_vals = 1:5,
    n_vals = seq(20L, 500L, 20L),
    n_reps = 30L,
    seed = NULL) {

  # Assertions
  assert_vector_cli(p_vals, min.len = 1, unique = TRUE)
  assert_numeric_cli(p_vals, lower = 1)
  assert_numeric_cli(n_vals, lower = 1)
  assert_vector_cli(n_vals, min.len = 1, unique = TRUE)
  assert_numeric_cli(n_reps, lower = 1)
  assert_numeric_cli(seed, lower = 0, null.ok = TRUE)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  cor_df  <- list()
  obs_df  <- list()
  cov_df  <- list()

  idx_cor <- 1
  idx_obs <- 1
  idx_cov <- 1

  for (p in p_vals) {
    for (n in n_vals) {

      ## Generate a fixed variance-covariance matrix for this (p,n) combination
      ## by simulating from a Wishart distribution
      # - The seed depends on p an n so that all replicates share the same
      #   Sigma matrix
      # - df = p + 2 ensures positive definiteness for var-cov matrices
      # - scale = diag(p) gives a baseline for variances
      set.seed(p * 10 + n)
      Sigma <- stats::rWishart(1, df = p + 2, Sigma = diag(p))[,,1]

      ## Simulate for each replicate
      for (rep in seq_len(n_reps)) {

        # Set the seed for this replicate
        set.seed(p * 10 + n * 10 + rep)

        # Generate data from MVN(vec(0), Sigma)
        smp <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)

        # Compute sample mean and covariance matrix
        mu_hat <- colMeans(smp)
        sigma_hat <- stats::cov(smp)

        # Compute Mahalanobis distance (squared) for each observation
        # based on the
        D2 <- stats::mahalanobis(smp, center = mu_hat, cov = sigma_hat)

        # Compute theoretical suitability (chi-squared)
        S_chisq <- 1 - stats::pchisq(D2, df = p)

        # Compute empirical suitability based on ECDF
        S_ecdf  <- 1 - stats::ecdf(D2)(D2)

        # Compute Pearson correlation coefficient between the two suitabilities
        # and store correlation result in main data frame
        cor_df[[idx_cor]] <- data.frame(
          p = p,
          n = n,
          rep = rep,
          cor = stats::cor(S_chisq, S_ecdf)
        )

        # Store correlation result in main data frame
        # Store individual observations
        obs_df[[idx_obs]] <- data.frame(
          p = p,
          n = n,
          rep = rep,
          D2 = D2,
          S_chisq = S_chisq,
          S_ecdf = S_ecdf
        )

        ## Compute variance-covariance matrix properties
        eig <- eigen(sigma_hat, symmetric = TRUE, only.values = TRUE)$values
        # Condition Number
        cond_num <- max(eig) / min(eig)  # measure of multicollinearity

        # Average absolute off-diagonal values (covariances/correlations)
        avg_abs_cor <- if (p > 1) {
          mean(abs(stats::cov2cor(sigma_hat)[lower.tri(sigma_hat)]))
        } else NA

        # Store properties
        cov_df[[idx_cov]] <- data.frame(
          p = p,
          n = n,
          cond_num = cond_num,
          avg_abs_cor = avg_abs_cor
        )

        idx_cor <- idx_cor + 1
        idx_obs <- idx_obs + 1
        idx_cov <- idx_cov + 1
      }
    }
  }

  # Bind results
  cor_df <- do.call(rbind, cor_df)
  obs_df <- do.call(rbind, obs_df)
  cov_df <- do.call(rbind, cov_df)

  # ------------------------------------------------------------------------------
  # Graph 1: Correlation between theoretical and empirical suitabilities
  # ------------------------------------------------------------------------------

  # Summarize correlations based on (n,p) combinations: mean and sd across
  # all replicates
  cor_summary <- cor_df |>
    dplyr::group_by(p, n) |>
    dplyr::summarise(
      mean_cor = mean(cor),
      sd_cor = stats::sd(cor),
      .groups = "drop"
    )

  # Plot: correlation vs. sample size, with error bars, by p
  cor_plot <- ggplot2::ggplot() +
    ggplot2::ylim(c(0.94, 1)) +
    ggplot2::geom_hline(yintercept = 1,
                        col = "black",
                        linetype = "dashed") +
    ggplot2::geom_point(data = cor_df,
                        ggplot2::aes(x = .data$n,
                                     y = .data$cor),
                        shape = 1,
                        col = "grey") +
    ggplot2::geom_pointrange(data = cor_summary,
                             ggplot2::aes(x = .data$n,
                                          y = .data$mean_cor,
                             ymin = .data$mean_cor - .data$sd_cor,
                             ymax = .data$mean_cor + .data$sd_cor),
                    size = 0.5,
                    linewidth = 0.5,
                    colour = "red") +
    ggplot2::facet_wrap(. ~ .data$p, nrow = 1, labeller = ggplot2::label_both) +
    ggplot2::labs(x = expression("Number of records (n)"),
                  y = "Pearson's correlation") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  # ------------------------------------------------------------------------------
  # Graph 2: Suitability vs. Mahalanobis distance
  # ------------------------------------------------------------------------------

  # Summarize observations across all replicates and sample sizes, then compute
  # the empirical suitability function on the pooled D2 values
  plot_df <- obs_df |>
    dplyr::group_by(p) |>
    dplyr::mutate(S_ecdf_pooled = 1 - stats::ecdf(D2)(D2)) |> # global ECDF per p
    dplyr::ungroup()

  # For each observation, we have both the S_chisq and the pooled S_ecdf
  # (use only a proportion of points)
  suit_plot <- ggplot2::ggplot(dplyr::sample_frac(plot_df, 0.1), ## sample 10% of points for plotting
                               ggplot2::aes(x = .data$D2,
                                            y = .data$S_ecdf_pooled)) +
    ggplot2::geom_point(alpha = 0.1,
                        ggplot2::aes(color = "ECDF"), size = 1.2) +
    ggplot2::stat_function(fun = function(x) 1 - stats::pchisq(x, df = 1),
                  aes(color = "Chi-squared")) +
    ggplot2::stat_function(fun = function(x) 1 - stats::pchisq(x, df = 2),
                  colour = "grey") +
    ggplot2::stat_function(fun = function(x) 1 - stats::pchisq(x, df = 3),
                  colour = "grey") +
    ggplot2::stat_function(fun = function(x) 1 - stats::pchisq(x, df = 4),
                  colour = "grey") +
    ggplot2::stat_function(fun = function(x) 1 - stats::pchisq(x, df = 5),
                  colour = "grey") +
    ggplot2::facet_wrap(. ~ .data$p, nrow = 1, labeller = ggplot2::label_both) +
    ggplot2::labs(x = expression("Squared Mahalanobis distance" ~ (D^2)),
                  y = "Habitat suitability") +
    ggplot2::scale_color_manual(name = NULL,
                       values = c("ECDF" = "red",
                                  "Chi-squared" = "grey")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
          panel.grid.minor = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(color = "black"))

  # ------------------------------------------------------------------------------
  # Graph 3: Mean correlation vs. condition number
  # ------------------------------------------------------------------------------

  # Merge the mean correlation with summary
  analysis_df <- merge(cor_summary, cov_df, by = c("p", "n"))

  # This measures the effect of collinearity on ECDF
  cond_plot <- ggplot2::ggplot(analysis_df,
                               ggplot2::aes(x = .data$cond_num,
                                            y = .data$mean_cor,
                                            colour = .data$n)) +
    ggplot2::ylim(c(0.94, 1)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(. ~ .data$p, nrow = 1, labeller = ggplot2::label_both) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Condition number (log)",
                  y = "Pearson's correlation",
                  color = "Num. records (n)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  # ------------------------------------------------------------------------------
  # Graph 4: Mean correlation vs. average absolute correlation (off-diagonal)
  # ------------------------------------------------------------------------------

  # Only for p > 1
  p_cor_avg <- ggplot2::ggplot(subset(analysis_df, analysis_df$p > 1),
                               ggplot2::aes(x = .data$avg_abs_cor,
                                            y = .data$mean_cor,
                                            color = .data$n)) +
    ggplot2::ylim(c(0.94, 1)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(. ~ .data$p, nrow = 1, labeller = ggplot2::label_both) +
    ggplot2::labs(x = "Average absolute off-diagonal correlation",
                  y = "Pearson's correlation",
                  color = "Num. records (n)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")


  # Return all results in a list
  l <- list(
    cor_plot = cor_plot,
    suit_plot = suit_plot,
    cond_plot = cond_plot,
    cor_avg_plot = p_cor_avg,
    cor_df = cor_df,
    obs_df = obs_df,
    cov_df = cov_df
  )

  return(l)
}
