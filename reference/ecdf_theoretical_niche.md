# Niche analysis using ECDF and chi-squared

Simulate niche suitability from Mahalanobis distance using both
chi-squared and empirical CDF transformations, for a given number of
predictor variables.

## Usage

``` r
ecdf_theoretical_niche(
  n,
  n_population = 10000L,
  sample_sizes = seq(20L, 500L, 20L),
  seed = NULL
)
```

## Arguments

- n:

  Integer; number of predictor variables (dimensions).

- n_population:

  Integer; size of simulated environmental population.

- sample_sizes:

  Integer vector of sample sizes to evaluate.

- seed:

  Optional integer seed for reproducibility.

## Value

A list with:

- corplot: ggplot object with correlation vs sample size.

- sample_data: matrix of simulated sample points.

- sample_niche: numeric vector of “true” niche suitability.

- chisq_suits: numeric vector, 1 - pchisq(Mahalanobis).

- ecdf_suits: numeric vector, 1 - ECDF(Mahalanobis).

- mahal_dists: numeric vector of Mahalanobis distances.

## Author

Luíz Fernando Esser

## Examples

``` r
# Create ECDF-niche based on personalized options:
n <- ecdf_theoretical_niche(n = 3,
                            n_population = 20000,
                            sample_sizes = seq(50, 1000, 50),
                            seed = 123)
```
