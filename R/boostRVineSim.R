#' Simulation from a Conditional R-Vine Copula Model
#'
#' This function simulates from an \code{boostRVineCopula} object.
#'
#' @param boostRVC an \code{boostRVineCopula} object to sample from.
#' @param U a matrix with d columns containing (randomly) uniformly distributed values in [0,1].
#' @param X matrix containing the covariates and having the same number of rows as \code{U}.
#'
#' @return A matrix with d columns containing the simulated data from the given d-dimensional conditional R-vine copula.
#'
#' @note The code of this function is mainly based on the R-package \code{\link{VineCopula}} by Nagler et al. but slightly adapted to the needs of this package.
#'
#' @examples
#' # load simulated data
#' data(data_vinecop)
#'
#' # create boostRVineMatrix
#' Matrix <- c(5, 0, 0, 0, 0,
#'             2, 2, 0, 0, 0,
#'             3, 3, 3, 0, 0,
#'             1, 4, 4, 4, 0,
#'             4, 1, 1, 1, 1)
#' Family <- c(0, 0, 0, 0, 0,
#'             404, 0, 0, 0, 0,
#'             1, 404, 0, 0, 0,
#'             401, 301, 304, 0, 0,
#'             1, 401, 304, 301, 0)
#' Matrix <- matrix(Matrix, 5, 5, byrow = TRUE)
#' Family <- matrix(Family, nrow = 5, ncol = 5, byrow = TRUE)
#' Formula <- matrix("~ .", 5, 5, byrow = TRUE)
#' boostRVM <- boostRVineMatrix(Matrix = Matrix,
#'                              family = Family,
#'                              formula = Formula)
#'
#' # fit object
#' object <- boostRVineSeqEst(U = data_vinecop[, 1:5],
#'                            X = data_vinecop[, -c(1:5)],
#'                            boostRVM = boostRVM,
#'                            control = boost_control(deselection = "attributable"),
#'                            cores = 1)
#' boostRVineSim(object,
#'               U = data_vinecop[, 1:5],
#'               X = data_vinecop[, -c(1:5)])
#'
#' @importFrom VineCopula BiCopHinv1 BiCopHinv2
#' @export
# The following codes are mainly based on the R-package VineCopula by Nagler et al.
boostRVineSim <- function(boostRVC, U, X) {

  if (nrow(U) != nrow(X)) {
    stop("Rows of U and X need to be of equal length!")
  }

  N <- nrow(U)
  o <- diag(boostRVC$boostRVM$Matrix)
  d <- length(o)
  boostRVC <- boostRVineNormalize(boostRVC)
  mat <- boostRVC$Matrix
  maxmat <- createMaxMat(mat)
  cindirect <- neededCondDistr(mat)$indirect

  revert <- function(m) {
    if (length(dim(m)) == 2) {
      return(m[nrow(m):1, ncol(m):1, drop = F])
    } else {
      return(m[nrow(m):1, ncol(m):1, , drop = F])
    }
  }
  revertlist <- function(l) {

    d <- length(l)
    out <- lapply(d:1, function(k) rev(l[[k]]))
    return(out)

  }
  pc <- revertlist(boostRVC$pair_copulas)
  mat <- revert(mat)
  maxmat <- revert(maxmat)
  cindirect <- revert(cindirect)

  vdirect <- vindirect <- array(dim = c(d, d, N))
  ## fill diagonal entries with independent uniforms
  id <- 1:d + d*(0:(d - 1)) + rep((0:(N - 1)) * d^2, each = d)
  vdirect[id] <- t(U[, rev(o)])
  vindirect[1, 1,] <- vdirect[1, 1,]

  for (i in 2:d) {
    for (k in (i - 1):1) {

      object <- pc[[k]][[i]]
      # parameter
      par1 <- predict(object, newdata = X, type = "parameter")$parameter
      par2 <- object$par2
      # family
      family <- getFams(as.numeric(object$family))
      fam <- rep(family[1], length(par1))
      if (length(family) != 1) {
        fam[par1 < 0] <- family[2]
      }

      m <- maxmat[k, i]
      u1 <- (if (mat[k, i] == m) vdirect else vindirect)[k, m,]
      vdirect[k, i,] <- BiCopHinv1(u1 = u1,
                                   u2 = vdirect[k + 1, i,],
                                   family = fam,
                                   par = par1,
                                   par2 = par2,
                                   check.pars = FALSE)
      if (i < d) {
        if (cindirect[k + 1, i]) {
          vindirect[k + 1, i,] <- BiCopHfunc2(u1 = u1,
                                              u2 = vdirect[k, i,],
                                              family = fam,
                                              par = par1,
                                              par2 = par2,
                                              check.pars = FALSE)
        }
      }

    }
  }

  out <- matrix(vdirect[1, , ], ncol = d, byrow = TRUE)
  # if (!is.null(boostRVC$Matrix)) {
  #   colnames(out) <- boostRVC$names
  # }
  out <- out[, sort(o[length(o):1], index.return = TRUE)$ix]

  return(out)

}


###########################################################
# Helper functions for boostRVineSim
###########################################################

boostRVineNormalize <- function(boostRVC) {

  oldOrder <- diag(boostRVC$boostRVM$Matrix)
  Matrix <- reorderRVineMatrix(boostRVC$boostRVM$Matrix)

  out <- list(Matrix = Matrix,
              pair_copulas = boostRVC$pair_copulas,
              names = rev(boostRVC$names[oldOrder]),
              covariates = boostRVC$covariates)

  return(out)

}
