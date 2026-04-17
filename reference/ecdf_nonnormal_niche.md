# Compare ECDF and Chi-square suitability under non-normal data

Script to run a simulation study to compare Chi-square vs. ECDF
approaches to quantify habitat suitability based on bivariate non-normal
data. Bivariate data was simulated based on environmental variables
(temperature and precipitation) using Gaussian copulas. Temperature
followed a normal distribution while precipitation followed a Weibull
distribution. The choices of the distributions were based on Haddad
(2021) - Theoretical and Applied Climatology (for temperature) and on
the estimation of rainfall in milimeters by Wilks (1989) - Journal of
Applied Meteorology. Because the relationship between temperature and
precipitation is complex across space (Rodrigo, 2022 - Theoretical and
Applied Climatology), we defined five correlation values between the two
variables.

## Usage

``` r
ecdf_nonnormal_niche(
  rho_vals = c(-0.7, -0.3, 0, 0.3, 0.7),
  n_vals = c(20L, 50L, 100L, 200L, 500L),
  n_reps = 10L,
  N_ref = 1e+05,
  shape_precip = NULL,
  scale_precip = NULL,
  mu_temp = NULL,
  sd_temp = NULL,
  seed = NULL
)
```

## Arguments

- rho_vals:

  Numeric vector; correlations between variables.

- n_vals:

  Integer vector; sample sizes.

- n_reps:

  Integer; number of replicates.

- N_ref:

  Integer; size of reference population for "true" parameters.

- shape_precip:

  Numeric; shape of the Weibull Distribution to obtain precipitation
  values.

- scale_precip:

  Numeric; scale of the Weibull Distribution to obtain precipitation
  values.

- mu_temp:

  Numeric; mean value to obtain normally distributed temperature values.

- sd_temp:

  Numeric; standard deviation to obtain normally distributed temperature
  values.

- seed:

  Optional integer for reproducibility.

## Value

A list with:

- suit_plot: ggplot of suitability vs Mahalanobis distance

- cor_df: correlation results

- obs_df: observation-level data

## Details

Simulates bivariate environmental data using Gaussian copulas with
non-normal marginals (Normal for temperature and Weibull for
precipitation), and evaluates agreement between chi-squared and ECDF
suitability.

## Author

Matheus T. Baumgartner

## Examples

``` r
# Create ECDF-niche based on personalized options:
n <- ecdf_nonnormal_niche(rho_vals = c(-0.7, -0.3, 0, 0.3, 0.7),
                          n_vals   = c(20L, 50L, 100L, 200L, 500L),
                          n_reps   = 10L,
                          N_ref    = 1e5,
                          seed     = 1991)
```
