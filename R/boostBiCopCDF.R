#' Distribution Function of a Conditional Bivariate Copula
#'
#' Cumulative distribution function (CDF) of an \code{boostBiCop} object.
#'
#' @param object an \code{boostBiCop} object.
#' @param U a matrix with two columns with values in [0,1].
#' @param X matrix containing the covariates and having the same number of rows as \code{U}.
#' @param ... unused.
#'
#' @return A numeric vector containing the distribution values of the considered \code{boostBiCop} object, evaluated at \code{U} and \code{X}.
#'
#' @importFrom VineCopula BiCopCDF
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
#' boostBiCopCDF(object,
#'               U = data_bicop[, 1:2],
#'               X = data_bicop[, -c(1:2)])
#'
#' @export

boostBiCopCDF <- function(object, U, X, ...) {

  if(missing(U) | missing(X)) {
    if (!"model_frame" %in% names(object)) {
      stop("Missing 'model_frame' in object due to 'light = TRUE'! Please provide arguments for 'U' and 'X'!")
    }
    U <- object$model_frame[, c(1, 2)]
    X <- object$model_frame[, -c(1, 2)]
  }

  if (ncol(U) != 2) {
    stop("U must have two columns!")
  }

  if (nrow(U) != nrow(X)) {
    stop("Rows of U and X do not have equal length!")
  }

  family <- as.numeric(object$family)
  par <- unlist(predict(object,
                        newdata = X,
                        type = "parameter"))
  par2 <- object$par2

  fam <- getFams(family)
  family <- rep(fam[1], length(par))
  if (length(fam) == 2) {
    family[par < 0] <- fam[2]
  }

  out <- BiCopCDF(u1 = U[, 1],
                  u2 = U[, 2],
                  family = family,
                  par = par,
                  par2 = par2,
                  check.pars = FALSE)

  return(out)

}
