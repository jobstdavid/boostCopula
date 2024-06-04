#' Generic functions
#'
#' @noRd
#' @export
print.boostBiCop <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  coef <- x$coefficients
  coef <- coef[coef != 0]

  cat(paste("Conditional bivariate copula: ", x$familyname, " (", x$family, ")", "\n", sep = ""))
  if (length(coef) != 0) {
    cat("Non-zero coefficients:\n")
    print.default(format(coef, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  }
  cat(paste("logLik: ", round(x$stats$loglik, 3),
            ", AIC: ", round(x$stats$aic, 3),
            ", degrees of freedom: ", x$stats$df,
            "\n", sep = ""))
  cat(paste("Iterations: ", x$iterations, sep = ""))

}
#' @noRd
#' @export
summary.boostBiCop <- function(object, digits = max(3, getOption("digits") - 3), ...) {

  print.boostBiCop(object, digits, ...)

}
#' @noRd
#' @export
print.boostRVineCopula <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat(paste0(ncol(x$boostRVM$Matrix),"-dimensional R-vine"))
  cat("\n")
  cat(paste("logLik: ", round(x$stats$loglik, 3),
            ", AIC: ", round(x$stats$aic, 3),
            ", degrees of freedom: ", x$stats$df,
            "\n", sep = ""))

}
#' @noRd
#' @export
summary.boostRVineCopula <- function(object, digits = max(3, getOption("digits") - 3), ...) {

  print.boostRVineCopula(object, digits, ...)

}
#' @noRd
#' @export
coef.boostBicop <- function(object, ...) {
  object$coefficients
}
#' @noRd
#' @export
logLik.boostBiCop <- function(object, ...) {

  structure(object$stats$loglik, df = object$stats$df, class = "logLik")

}
#' @noRd
#' @export
AIC.boostBiCop <- function(object, ...) {
  object$stats$aic
}
#' @noRd
#' @export
logLik.boostRVineCopula <- function(object, ...) {

  structure(object$stats$loglik, df = object$stats$df, class = "logLik")

}
#' @noRd
#' @export
AIC.boostRVineCopula <- function(object, ...) {
  object$stats$aic
}
