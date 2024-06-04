#' Control Parameters for the Gradient-Boosting Algorithm
#'
#' @param maxit integer; giving the number of initial boosting iterations. Default: \code{maxit = 500}.
#' @param nu double; defining the step size or shrinkage parameter. Default: \code{nu = 0.1}.
#' @param mstop character; stopping criterion "\code{max}", "\code{aic}" (default), "\code{cv}".
#' @param nfolds integer; number of folds in cross validation \code{mstop = "cv"}.
#' @param foldid vector; an optional vector of values between 1 and \code{nfolds} identifying the fold for each observation.
#' Default: \code{foldid = NULL}.
#' @param deselection character; "\code{attributable}", "\code{cumulative}" or "\code{none}" (default).
#' @param gamma double; defining the threshold for \code{deselection}. Default: \code{gamma = 0.01}.
#' @param center logical; indicating if the covariates should be mean centered \code{TRUE} (default) before fitting or not \code{FALSE}.
#' @param cores.control integer; the number of cores used for computations. Default: \code{cores = 1}.
#'
#' @return List with specified parameters.
#'
#' @export
boost_control <- function(maxit = 500,
                          nu = 0.1,
                          mstop = "aic",
                          nfolds = 10,
                          foldid = NULL,
                          deselection = "none",
                          gamma = 0.01,
                          center = TRUE,
                          cores.control = 1) {

  out <- list(maxit = maxit,
              nu = nu,
              mstop = mstop,
              nfolds = nfolds,
              foldid = foldid,
              deselection = deselection,
              gamma = gamma,
              center = center,
              cores.control = cores.control)

  class(out) <- "boost_control"
  return(out)

}

