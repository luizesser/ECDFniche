# Create distance–suitability plot

Create distance–suitability plot

## Usage

``` r
create_distance_suitability_plot(analysis_results)
```

## Arguments

- analysis_results:

  List returned by
  [`ecdf_theoretical_niche()`](https://luizesser.github.io/ECDFniche/reference/ecdf_theoretical_niche.md).

## Value

A ggplot object.

## Examples

``` r
# Create ECDF-niche based on personalized options:
res <- ecdf_theoretical_niche(n = 3,
                              n_population = 20000,
                              sample_sizes = seq(50, 1000, 50),
                              seed = 123)

# Plot analysis results
create_distance_suitability_plot(res)

```
