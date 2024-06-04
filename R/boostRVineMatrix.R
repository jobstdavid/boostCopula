#' Conditional R-Vine Copula Matrix Notation
#'
#' This function creates an \code{boostRVineMatrix} object for a d-dimensional conditional R-vine copula
#' containing the R-vine structure, the specified pair-copula families and the corresponding formulas.
#'
#' @param Matrix Lower triangular d x d integer matrix that defines the R-vine tree structure.
#' @param family Lower triangular d x d integer matrix that assigns the pair-copula families to each (conditional) pair defined by Matrix
#' (default: \code{family = array(NA, dim = dim(Matrix))}). The available bivariate copula families are given in \code{\link{boostBiCopEst}}.
#' @param formula Lower triangular d x d character matrix that assigns the formula to the corresponding copula family specified in \code{family}
#' (default: \code{formula = array(NA, dim = dim(Matrix))}).
#'
#' @return An \code{boostRVineMatrix} object.
#'
#' @export
boostRVineMatrix <- function(Matrix,
                             family = array(NA, dim = dim(Matrix)),
                             formula = array(NA, dim = dim(Matrix))) {

  out <- list(Matrix = Matrix,
              family = family,
              formula = formula)
  class(out) <- "boostRVineMatrix"
  return(out)

}








