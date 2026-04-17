#' Mahalanobis Distance Classifier for Ecological Niche Modeling
#'
#' A custom \code{caret} model specification implementing a Mahalanobis
#' distance-based classifier for ecological niche modeling (ENM) and
#' species distribution modeling (SDM). This implementation supports both
#' parametric (chi-squared) and nonparametric (empirical cumulative
#' distribution function; ECDF) transformations of Mahalanobis distances
#' into suitability scores.
#'
#' The model is trained using presence-only data to estimate the centroid
#' and covariance structure of environmental conditions associated with
#' species occurrences. Suitability is then derived as the inverse tail
#' probability of the Mahalanobis distance between new observations and
#' the estimated niche centroid.
#'
#' Two approaches are available to transform Mahalanobis distances into
#' probabilities:
#' \itemize{
#'   \item \code{"chisq"}: assumes distances follow a chi-squared
#'   distribution with degrees of freedom equal to the number of predictors.
#'   \item \code{"ecdf"}: uses the empirical cumulative distribution
#'   function of training distances, providing a nonparametric estimate
#'   of suitability.
#' }
#'
#' The ECDF-based approach is particularly useful when the assumption of
#' multivariate normality is violated, which is common in ecological data.
#'
#' This model can be used within the \code{caret::train()} framework,
#' enabling resampling, tuning, and ensemble modeling workflows for
#' ecological niche modeling.
#'
#' @section Model Parameters:
#' \describe{
#'   \item{abs}{Logical. If \code{TRUE}, predictions are binarized using a
#'   fixed threshold (default: 0.05). If \code{FALSE}, the class with the
#'   highest predicted probability is returned.}
#'
#'   \item{method}{Character. Method used to convert Mahalanobis distances
#'   into suitability values. Options are \code{"chisq"} or \code{"ecdf"}.}
#' }
#'
#' @section Details:
#' The Mahalanobis distance defines an ellipsoidal niche in environmental
#' space. Under the chi-squared formulation, suitability decreases as the
#' distance from the niche centroid increases. The ECDF formulation
#' relaxes distributional assumptions by estimating suitability directly
#' from the empirical distribution of distances observed in presence data.
#'
#' Predictions return class probabilities for \code{"presence"} and
#' \code{"pseudoabsence"}, allowing flexible thresholding and ensemble
#' integration.
#'
#' @section Usage in caret:
#' This object can be supplied to \code{caret::train()} as a custom model:
#'
#' \preformatted{
#' library(caret)
#'
#' model <- train(
#'   x = predictors,
#'   y = response,
#'   method = mahal.dist,
#'   trControl = trainControl(classProbs = TRUE)
#' )
#' }
#'
#' You can also run only ECDF by adjusting the tuning grid:
#' \preformatted{
#' library(caret)
#'
#' grid <- expand.grid(
#'   abs = c(TRUE, FALSE),
#'   method = "ecdf"
#' )
#'
#' model <- train(
#'   x = predictors,
#'   y = response,
#'   method = mahal.dist,
#'   tuneGrid = grid,
#'   trControl = trainControl(classProbs = TRUE)
#' )
#' }
#'
#' @seealso \code{\link[stats]{mahalanobis}}, \code{\link[stats]{ecdf}},
#'   \code{\link[caret]{train}}
#'
#' @keywords models niche-modeling species-distribution-modeling mahalanobis
#'
#' @importFrom stats cov mahalanobis pchisq ecdf
#' @export
mahal.dist <- list(

  label = "Mahalanobis Distance Classifier",
  library = NULL,
  type = "Classification",

  parameters = data.frame(
    parameter = c("abs", "method"),
    class = c("logical", "character"),
    label = c("Absolute Binarization", "Suitability method")
  ),

  grid = function(x, y, len = NULL, search = "grid") {
    expand.grid(
      abs = c(TRUE, FALSE),
      method = c("chisq", "ecdf")
    )
  },

  fit = function(x, y, wts, param, lev, last, classProbs, ...) {

    # Use only presence data
    presence_data <- x[y == "presence", , drop = FALSE]

    if (nrow(presence_data) < 2) {
      stop("Not enough 'presence' data points to calculate covariance.")
    }

    # Core Mahalanobis parameters
    center_vec <- colMeans(presence_data, na.rm = TRUE)
    cov_mat <- stats::cov(presence_data)
    inv_cov_matrix <- solve(cov_mat)

    # Training distances (for ECDF)
    d2_train <- stats::mahalanobis(
      x = presence_data,
      center = center_vec,
      cov = inv_cov_matrix,
      inverted = TRUE
    )

    # Store model
    result <- list(
      center = center_vec,
      inv_cov = inv_cov_matrix,
      df = ncol(x),
      abs = param$abs,
      method = param$method,
      ecdf_fun = stats::ecdf(d2_train),
      levels = lev
    )

    return(result)
  },

  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {

    probs <- .mahal.dist$prob(modelFit, newdata)

    if (modelFit$abs) {
      pred <- ifelse(
        probs[, modelFit$levels[1]] >= 0.05,
        modelFit$levels[1],
        modelFit$levels[2]
      )
    } else {
      pred <- colnames(probs)[apply(probs, 1, which.max)]
    }

    factor(pred, levels = modelFit$levels)
  },

  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {

    d2 <- stats::mahalanobis(
      x = newdata,
      center = modelFit$center,
      cov = modelFit$inv_cov,
      inverted = TRUE
    )

    # Switch between methods
    if (modelFit$method == "chisq") {

      p_presence <- 1 - stats::pchisq(d2, df = modelFit$df)

    } else if (modelFit$method == "ecdf") {

      p_presence <- 1 - modelFit$ecdf_fun(d2)

    } else {
      stop("Unknown method: ", modelFit$method)
    }

    prob_df <- data.frame(
      presence = p_presence,
      pseudoabsence = 1 - p_presence
    )

    colnames(prob_df) <- modelFit$levels
    return(prob_df)
  },

  predictors = function(x, ...) {
    names(x$center)
  },

  varImp = function(object, ...) {
    # simple proxy: variance contribution
    data.frame(Overall = diag(solve(object$inv_cov)))
  },

  levels = function(x) x$levels,

  tags = c("mahalanobis", "Presence-Only", "ECDF", "Nonparametric")
)
