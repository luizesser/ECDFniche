# Run full ECDF–Mahalanobis analysis

Convenience function that reproduces the three figures from the original
manuscript for 1–5 dimensions.

## Usage

``` r
run_ecdf_mahal_analysis(dims = 1:5, seed = 3L)
```

## Arguments

- dims:

  Integer vector of dimensions (default 1:5).

- seed:

  Optional seed for reproducibility.

## Value

A list containing:

- analyses: list of ecdf_theoretical_niche() outputs.

- figure1, figure2, figure3: grobs with arranged plots.

## Examples

``` r
# Recreate original manuscript output:
set.seed(3)
full_res <- run_ecdf_mahal_analysis(dims = 1:5)



```
