#' Compare ECDF and Chi-square suitability under non-normal data
#'
#' Simulates bivariate environmental data using Gaussian copulas with
#' non-normal marginals (Normal for temperature and Weibull for precipitation),
#' and evaluates agreement between chi-squared and ECDF suitability.
#'
#' @param rho_vals Numeric vector; correlations between variables.
#' @param n_vals Integer vector; sample sizes.
#' @param n_reps Integer; number of replicates.
#' @param N_ref Integer; size of reference population for "true" parameters.
#' @param seed Optional integer for reproducibility.
#' @param mu_temp  Numeric; mean value to obtain normally distributed temperature values.
#' @param sd_temp  Numeric; standard deviation to obtain normally distributed temperature values.
#' @param shape_precip Numeric; shape of the Weibull Distribution to obtain precipitation values.
#' @param scale_precip  Numeric; scale of the Weibull Distribution to obtain precipitation values.
#'
#' @description
#' Script to run a simulation study to compare Chi-square vs. ECDF approaches
#' to quantify habitat suitability based on bivariate non-normal data.
#' Bivariate data was simulated based on environmental variables (temperature
#' and precipitation) using Gaussian copulas. Temperature followed a normal
#' distribution while precipitation followed a Weibull distribution.
#' The choices of the distributions were based on Haddad (2021) - Theoretical
#' and Applied Climatology (for temperature) and on the estimation of rainfall
#' in milimeters by Wilks (1989) - Journal of Applied Meteorology.
#' Because the relationship between temperature and precipitation is complex
#' across space (Rodrigo, 2022 - Theoretical and Applied Climatology), we
#' defined five correlation values between the two variables.
#'
#' @return A list with:
#' \itemize{
#'   \item suit_plot: ggplot of suitability vs Mahalanobis distance
#'   \item cor_df: correlation results
#'   \item obs_df: observation-level data
#' }
#'
#' @author Matheus T. Baumgartner
#'
#' @examples
#' # Create ECDF-niche based on personalized options:
#' n <- ecdf_nonnormal_niche(rho_vals = c(-0.7, -0.3, 0, 0.3, 0.7),
#'                           n_vals   = c(20L, 50L, 100L, 200L, 500L),
#'                           n_reps   = 10L,
#'                           N_ref    = 1e5,
#'                           seed     = 1991)
#'
#' @importFrom MASS mvrnorm
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual geom_vline geom_histogram
#'             theme_void theme element_rect labs facet_wrap theme_bw label_both xlim unit
#' @importFrom ggpp geom_plot_npc
#' @importFrom dplyr mutate group_by ungroup
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @importFrom stats pnorm qnorm qweibull cov cov2cor mahalanobis pchisq ecdf cor
#' @import checkCLI
#'
#' @global .data data
#'
#' @export
ecdf_nonnormal_niche <- function(
    rho_vals = c(-0.7, -0.3, 0, 0.3, 0.7),
    n_vals   = c(20L, 50L, 100L, 200L, 500L),
    n_reps   = 10L,
    N_ref    = 1e5,
    shape_precip = NULL,
    scale_precip = NULL,
    mu_temp = NULL,
    sd_temp = NULL,
    seed     = NULL) {

  # Assertions
  assert_vector_cli(rho_vals, min.len = 1, unique = TRUE)
  assert_numeric_cli(rho_vals, lower = -1, upper = 1)
  assert_numeric_cli(n_vals, lower = 1)
  assert_vector_cli(n_vals, min.len = 1, unique = TRUE)
  assert_numeric_cli(n_reps, lower = 1)
  assert_numeric_cli(N_ref, lower = 1L, null.ok = TRUE)
  assert_numeric_cli(seed, lower = 0, null.ok = TRUE)

  assert_numeric_cli(mu_temp, len = 1, null.ok = TRUE)
  assert_numeric_cli(sd_temp, len = 1, null.ok = TRUE)
  assert_numeric_cli(shape_precip, len = 1, null.ok = TRUE)
  assert_numeric_cli(scale_precip, len = 1, null.ok = TRUE)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ------------------------------------------------------------------------------
  # 1. Set up the true distribution parameters (fixed across simulations)
  # ------------------------------------------------------------------------------
  ## Marginals
  # Temperature ~ Normal(mean = 20, sd = 5)
  mu_temp <- ifelse(is.null(mu_temp), 20, mu_temp)
  sd_temp <- ifelse(is.null(sd_temp), 5,  sd_temp)

  # Precipitation ~ Weibull(shape = 2, scale = 10)
  shape_precip <- ifelse(is.null(shape_precip), 2,  shape_precip)
  scale_precip <- ifelse(is.null(scale_precip), 10, scale_precip)

  # ------------------------------------------------------------------------------
  # 2. Generate a huge reference sample to obtain the 'true' population parameters
  #    (mean vector and covariance matrix) for each value of rho
  # ------------------------------------------------------------------------------

  # Store 'true' parameters in a list indexed by rho values
  true_params <- list()


  for (rho in rho_vals) {
    ## Generate reference data
    # Specify a variance-covariance matrix for the Gaussian copula
    Sigma_cop <- matrix(c(1, rho, rho, 1), 2)

    # Generate bivariate normal data for each variable based on specified correlation
    z_ref <- MASS::mvrnorm(N_ref, mu = c(0, 0), Sigma = Sigma_cop)

    # Transform data into Cumulative Distribution: Uniform(0, 1)
    u_ref <- stats::pnorm(z_ref)

    # Transform the copula into marginal distributions
    temp_ref   <- stats::qnorm(u_ref[,1], mean = mu_temp, sd = sd_temp)
    precip_ref <- stats::qweibull(u_ref[,2], shape = shape_precip, scale = scale_precip)

    X_ref <- cbind(temp_ref, precip_ref)

    # Calculate 'true' correlation matrix
    cov_true <- stats::cov(X_ref)
    cor_true <- stats::cov2cor(cov_true)

    true_params[[as.character(rho)]] <- list(
      mu    = colMeans(X_ref), # Calculate 'true' mean vector using marginal means
      Sigma = stats::cov(X_ref),  # Calculate 'true' variance-covariance matrix
      Cor   = cor_true           # Store 'true' correlation matrix
    )
  }


  # ------------------------------------------------------------------------------
  # 3. Simulation loop over rho, n, and replicates
  # ------------------------------------------------------------------------------

  # List to store correlation results
  cor_list <- list()
  obs_list <- list()
  rep_data <- list()

  idx_cor <- 1
  idx_obs <- 1

  for (rho in rho_vals) {
    # Generate sample from the copula model with the same rho
    Sigma_cop <- matrix(c(1, rho, rho, 1), 2)

    # Retrieve true parameters for this rho
    mu_true  <- true_params[[as.character(rho)]]$mu
    cov_true <- true_params[[as.character(rho)]]$Sigma

    for (n in n_vals) {
      for (rep in seq_len(n_reps)) {

        z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma_cop) # generate bivariate normal data for copula
        u <- stats::pnorm(z) # convert bivariate normal into cumulative distirbution

        temp   <- stats::qnorm(u[,1], mean = mu_temp, sd = sd_temp) # generate marginal temperature
        precip <- stats::qweibull(u[,2], shape = shape_precip, scale = scale_precip) # generate marginal precipitation

        X <- cbind(temp, precip) # combine into data frame

        # Calculate squared Mahalanobis distances using 'true' population parameters
        # Using 'true' values instead of sample values isolates the effect of non-normality
        # from estimation error
        D2 <- stats::mahalanobis(X, center = mu_true, cov = cov_true)

        # Calculate 'theoretical' suitability (based on chi-square)
        # This is 'wrong' here because distribution is non-normal, but to keep as a comparison
        S_chisq <- 1 - stats::pchisq(D2, df = 2)

        # Calculate 'empirical' suitability (based on ECDF)
        S_ecdf  <- 1 - stats::ecdf(D2)(D2)

        cor_list[[idx_cor]] <- data.frame(
          rho = rho,
          n   = n,
          rep = rep,
          cor = stats::cor(S_chisq, S_ecdf) # Calculate Pearson's correlation between the two suitabilities
        )

        obs_list[[idx_obs]] <- data.frame(
          rho = rho,
          n   = n,
          rep = rep,
          D2  = D2,
          S_chisq = S_chisq,
          S_ecdf  = S_ecdf,
          diff = S_ecdf - S_chisq
        )

        idx_cor <- idx_cor + 1
        idx_obs <- idx_obs + 1

        # Save the replicate for plotting
        rep_data[[paste0("rho", rho, "_n", n, "_rep", rep)]] <- data.frame(
          rho = rho,
          n = n,
          rep = rep,
          D2 = D2,
          S_chisq = S_chisq,
          S_ecdf = S_ecdf
          )
      }
    }
  }

  # ------------------------------------------------------------------------------
  # 4. Build graph
  # ------------------------------------------------------------------------------

  # Combine all results in a single data frame
  cor_df <- do.call(rbind, cor_list)
  obs_df <- do.call(rbind, obs_list)
  complete_obs <- do.call(rbind, rep_data)

  # Calculate difference between ECDF and Chi-square
  complete_obs <- complete_obs |>
    dplyr::mutate(diff = .data$S_ecdf - .data$S_chisq)

  # Create a data frame for inset histograms
  inset_plots <- complete_obs |>
    dplyr::group_by(n) |>
    tidyr::nest() |>
    dplyr::mutate(plot = purrr::map(data, ~ ggplot2::ggplot(.x, ggplot2::aes(x = .data$diff)) +
                                      ggplot2::xlim(c(-0.3, 0.3)) +
                                      ggplot2::geom_histogram(bins = 20, fill = "gray80", color = NA) +
                                      ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
                                      ggplot2::theme_void() + # removes axes/background for a clean inset
                                      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white",
                                                                                             color = "black",
                                                                                             linewidth = 0.3))))

  # Create final plot
  suit_plot <- ggplot2::ggplot(complete_obs) +
    # Print suitability vs. Mahalanobis distance
    ggplot2::geom_point(ggplot2::aes(x = D2, y = S_ecdf, color = factor(rho)), size = 1, shape = 21) +
    ggplot2::scale_color_viridis_d(direction = -1) +
    ggplot2::geom_point(ggplot2::aes(x = .data$D2, y = .data$S_chisq), color = "red", size = 1, shape = 21) +

    # Print inset histograms
    ggpp::geom_plot_npc(data = inset_plots,
                        ggplot2::aes(npcx = 0.95, npcy = 0.65, label = plot),
                        vp.width = 0.3, vp.height = 0.25,
                        hjust = 1, vjust = 0) +

    ggplot2::facet_wrap(. ~ n, nrow = 1, labeller = ggplot2::label_both) +
    ggplot2::labs(color = expression(rho),
         x = expression("Squared Mahalanobis distance" ~ (D^2)),
         y = "Habitat Suitability") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.key.size = ggplot2::unit(3, "point"))

  l <- list(
    suit_plot = suit_plot,
    cor_df = cor_df,
    obs_df = obs_df
  )
  return(l)
}
