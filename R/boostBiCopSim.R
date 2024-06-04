#' Simulation from a Conditional Bivariate Copula
#'
#' Simulation from a \code{boostBiCop} object.
#'
#' @param object an \code{boostBiCop} object.
#' @param N integer, specifying the number of simulation samples.
#' @param X matrix containing the covariates.
#' @param ... unused.
#'
#' @return A matrix with two columns containing the samples of the \code{boostBiCop} object, evaluated at \code{X}.
#'
#' @importFrom VineCopula BiCopSim
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
#' boostBiCopSim(object,
#'               N = 10,
#'               X = data_bicop[, -c(1:2)])
#'
#' @export
boostBiCopSim <- function(object, N, X, ...) {

  if (missing(X)) {
    if (!"model_frame" %in% names(object)) {
      stop("Missing 'model_frame' in object due to 'light = TRUE'! Please provide an argument for 'X'!")
    }
    X <- object$model_frame[, -c(1, 2)]
  }
  N_tmp <- nrow(X)
  if ((N < N_tmp) || (N > N_tmp)) {
    X <- X[sample(N_tmp, N, replace = TRUE), ]
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

  out <- BiCopSim(N = N,
                  family = family,
                  par = par,
                  par2 = par2,
                  check.pars = FALSE)

  return(out)

}
