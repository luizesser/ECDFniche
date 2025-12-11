# ECDFniche

<!-- badges: start -->
<!-- badges: end -->

Luíz Fernando Esser

ECDFniche is an R package for converting environmental distances into probability-like values using empirical cumulative distribution functions (ECDF), aimed at ecological niche and species distribution modeling workflows.

## Package overview

ECDFniche implements an ECDF-based mapping from distances in environmental space to values on a 0–1 scale that behave like cumulative probabilities.
The package is designed primarily for presence-only niche modeling, where distance-based algorithms (e.g., Mahalanobis, Euclidean or other metrics) need a consistent way to convert distances into suitability or probability-like surfaces.
The package is motivated by the classic use of Mahalanobis distance and its mapping to a chi-squared distribution to obtain probabilities, while removing the need for strict parametric assumptions.
By relying on the empirical cumulative distribution of observed distances, ECDFniche can be used with multiple distance metrics and remains flexible when the theoretical distribution of the distance is unknown or violated.

## Installation

You can install the development version of ECDFniche from [GitHub](https://github.com/luizesser/ECDFniche) with:

``` r
install.packages("devtools")
devtools::install_github("luizesser/caretSDM")
```

The package is also available on CRAN. Users are able to install it using the following code:

``` r
install.packages("caretSDM")
```

# Citation

If you use ECDFniche in a publication, please cite the associated methodological paper on ECDF-based distance-to-suitability mapping in ecological niche modeling once available.
