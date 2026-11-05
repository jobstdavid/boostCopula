#' Predict from Conditional Bivariate Copula
#'
#' @param object an object of class \code{boostBiCop}.
#' @param newdata data frame of covariates for which to predict.
#' @param type character; either "\code{parameter}" (default) for the copula parameter,
#'"\code{tau}" for Kendall's Tau or "\code{response}" for the linear predictor.
#' @param ... unused.
#'
#' @return a data frame.
#'
#' @examples
#' # load simulated data
#' data(data_bicop)
#'
#' # fit object
#' object <- boostBiCopEst(formula = ~.,
#'                         U = data_bicop[, 1:2],
#'                         X = data_bicop[, -c(1:2)],
#'                         family = 301,
#'                         control = boost_control(deselection = "attributable"))
#' predict(object, type = "parameter")
#' predict(object, newdata = data_bicop, type = "tau")
#'
#' @importFrom VineCopula BiCopTau2Par
#' @export
predict.boostBiCop <- function(object, newdata, type = "parameter", ...) {

  if (!type %in% c("parameter", "tau", "response")) {
    stop("'type' should be one of 'parameter', 'tau' or 'response'!")
  }

  # subset data frame
  if (!missing(newdata)) {

    covs <- names(object$coefficients)
    # `setdiff`, not `covs[-index]`: a model without an intercept has no index
    # to drop, and `-integer(0)` would select no covariate at all.
    vars <- setdiff(covs, "(Intercept)")

    # check model matrix (faster than model.matrix)
    if (any(!vars %in% colnames(newdata))) {
        stop("Covariates in 'newdata' do not match with covariates in 'coefficients'!")
    }
    X <- as.matrix(newdata[, vars, drop = FALSE])

    if ("(Intercept)" %in% covs) {
      X <- cbind("(Intercept)" = 1, X)
    }
    X <- X[, covs, drop = FALSE]

  } else {
    if (!"model_frame" %in% names(object)) {
      stop("Missing 'model_frame' in object due to 'light = TRUE'! Please provide an argument for 'newdata'!")
    }
    X <- as.matrix(object$model_frame[, -c(1, 2)])
  }

  family <- as.numeric(object$family)
  coef <- object$coefficients
  # generalized linear model
  out <- matvecmult_eigen(X, coef, 1)

  if (type == "tau" | type == "parameter") {

    # Kendall's tau
    out <- tanh(out)

    # copula parameter
    if (type == "parameter") {
      out <- tau2par(family, out)
    }

  }

  out <- data.frame(out)
  colnames(out) <- type

  return(out)

}



