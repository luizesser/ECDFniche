# Simulations and analyses of Mahalanobis distance-based habitat suitability

The objective is to compare the performance of habitat suitability
calculated based on chi-squared cumulative distribution function and
Empirical Cumulative Distribution Function (ECDF)

## Usage

``` r
ecdf_compare_niche(
  p_vals = 1:5,
  n_vals = seq(20L, 500L, 20L),
  n_reps = 30L,
  seed = NULL
)
```

## Arguments

- p_vals:

  Integer vector; number of predictor variables (dimensions).

- n_vals:

  Integer vector; number of records (sample sizes).

- n_reps:

  Integer; number of replicates per combination.

- seed:

  Optional integer for reproducibility.

## Value

A list with:

- cor_plot: ggplot of correlation vs sample size.

- suit_plot: ggplot of suitability vs Mahalanobis distance.

- cond_plot: ggplot of correlation vs condition number.

- cor_df: raw correlation data.

- obs_df: observation-level data.

- cov_df: covariance diagnostics.

## Details

Performs replicated simulations of multivariate normal data to evaluate
the agreement between suitability derived from chi-squared distribution
and empirical cumulative distribution function (ECDF).

## Author

Matheus T. Baumgartner

## Examples

``` r
# Create ECDF-niche based on personalized options:
n <- ecdf_compare_niche(p_vals = 1:3,
                        n_vals = seq(50L, 500L, 50L),
                        n_reps = 10L,
                        seed = 1991)
```
