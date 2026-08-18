# Mahalanobis Distance Classifier for Ecological Niche Modeling

A custom `caret` model specification implementing a Mahalanobis
distance-based classifier for ecological niche modeling (ENM) and
species distribution modeling (SDM). This implementation supports both
parametric (chi-squared) and nonparametric (empirical cumulative
distribution function; ECDF) transformations of Mahalanobis distances
into suitability scores.

## Usage

``` r
mahal.dist
```

## Details

The model is trained using presence-only data to estimate the centroid
and covariance structure of environmental conditions associated with
species occurrences. Suitability is then derived as the inverse tail
probability of the Mahalanobis distance between new observations and the
estimated niche centroid.

Two approaches are available to transform Mahalanobis distances into
probabilities:

- `"chisq"`: assumes distances follow a chi-squared distribution with
  degrees of freedom equal to the number of predictors.

- `"ecdf"`: uses the empirical cumulative distribution function of
  training distances, providing a nonparametric estimate of suitability.

The ECDF-based approach is particularly useful when the assumption of
multivariate normality is violated, which is common in ecological data.

This model can be used within the `caret::train()` framework, enabling
resampling, tuning, and ensemble modeling workflows for ecological niche
modeling.

## Model Parameters

- abs:

  Logical. If `TRUE`, predictions are binarized using a fixed threshold
  (default: 0.05). If `FALSE`, the class with the highest predicted
  probability is returned.

- method:

  Character. Method used to convert Mahalanobis distances into
  suitability values. Options are `"chisq"` or `"ecdf"`.

## Details

The Mahalanobis distance defines an ellipsoidal niche in environmental
space. Under the chi-squared formulation, suitability decreases as the
distance from the niche centroid increases. The ECDF formulation relaxes
distributional assumptions by estimating suitability directly from the
empirical distribution of distances observed in presence data.

Predictions return class probabilities for `"presence"` and
`"pseudoabsence"`, allowing flexible thresholding and ensemble
integration.

## Usage in caret

This object can be supplied to `caret::train()` as a custom model:


    library(caret)

    model <- train(
      x = predictors,
      y = response,
      method = mahal.dist,
      trControl = trainControl(classProbs = TRUE)
    )

You can also run only ECDF by adjusting the tuning grid:


    library(caret)

    grid <- expand.grid(
      abs = c(TRUE, FALSE),
      method = "ecdf"
    )

    model <- train(
      x = predictors,
      y = response,
      method = mahal.dist,
      tuneGrid = grid,
      trControl = trainControl(classProbs = TRUE)
    )

## See also

[`mahalanobis`](https://rdrr.io/r/stats/mahalanobis.html),
[`ecdf`](https://rdrr.io/r/stats/ecdf.html), `train`
