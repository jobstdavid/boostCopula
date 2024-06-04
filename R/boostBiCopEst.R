#' Gradient-Boosted Estimation of GLMs for Conditional Bivariate Copulas
#'
#' Gradient-boosted estimation of GLMs for conditional bivariate copulas.
#'
#' @param formula  an formula object, starting with the
#' \code{~} operator and followed by the covariates, which are separated by \code{+} operators.
#' @param U a matrix with two columns with values in [0,1].
#' @param X matrix containing the covariates and having the same number of rows as \code{U}.
#' @param family integer; defining the bivariate copula family:\cr
#' 1 = Gaussian copula \cr
#' 301 = double Clayton I copula (0 and 90 degrees) \cr
#' 302 = double Clayton II copula (0 and 270 degrees) \cr
#' 303 = double Clayton III copula (180 and 90 degrees) \cr
#' 304 = double Clayton IV copula (180 and 270 degrees) \cr
#' 401 = double Gumbel I copula (0 and 90 degrees) \cr
#' 402 = double Gumbel II copula (0 and 270 degrees) \cr
#' 403 = double Gumbel III copula (180 and 90 degrees) \cr
#' 404 = double Gumbel IV copula (180 and 270 degrees) \cr
#' @param control control parameters for the gradient-boosting estimation.
#' The \code{\link{boost_control}} function provides the default setting.
#' @param na.action a function which indicates what should happen if the data contains \code{NA}s. Default: \code{na.omit}.
#' @param light logical; if \code{FALSE} (default) the model frame is kept in the \code{boostBiCop} object.
#' Otherwise, if \code{TRUE}, the model frame is omitted.
#' @param ... unused.
#'
#' @return An \code{boostBiCop} object.
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
#' summary(object)
#' plot(object)
#'
#' @export
boostBiCopEst <- function(formula,
                          U,
                          X,
                          family,
                          control = boost_control(),
                          na.action = na.omit,
                          light = FALSE,
                          ...) {

  object <- boostBiCopSelect(formula = formula,
                             U = U,
                             X = X,
                             familyset = family,
                             selectioncrit = "loglik",
                             indeptest = NA,
                             control = control,
                             na.action = na.action,
                             light = light,
                             cores = 1,
                             ...)

  return(object)


}


