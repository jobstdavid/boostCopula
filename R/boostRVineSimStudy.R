#' Simulation from an Artificial Conditional R-Vine Copula Model
#'
#' This function simulates from a given d-dimensional artificial conditional R-vine copula model.
#'
#' @param boostRVC an \code{boostRVineCopula} object to sample from.
#' @param U a matrix with d columns containing (randomly) uniformly distributed values in [0,1].
#' @param response vector containing the linear predictor assumed among all bivariate conditional copulas.
#'
#' @return A matrix with d columns containing the simulated data from the given d-dimensional conditional R-vine copula.
#'
#' @note The code of this function is mainly based on the R-package \code{\link{VineCopula}} by Nagler et al. but slightly adapted to the needs of this package.
#'
#' @export
boostRVineSimStudy <- function(boostRVC, U, response) {

  if (nrow(U) != length(response)) {
    stop("Rows of 'U' and length of 'response' need to be equal!")
  }

  N <- nrow(U)
  o <- diag(boostRVC$Matrix)
  d <- length(o)
  boostRVC <- Normalize(boostRVC)
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
  family <- revert(boostRVC$family)
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

      # parameter
      par1 <- tau2par(family[k, i], tanh(response))
      par2 <- 0
      # family
      fams <- getFams(as.numeric(family[k, i]))
      fam <- rep(fams[1], length(par1))
      if (length(fams) != 1) {
        fam[par1 < 0] <- fams[2]
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
  if (!is.null(boostRVC$names)) {
    colnames(out) <- boostRVC$names
  }
  out <- out[, sort(o[length(o):1], index.return = TRUE)$ix]

  return(out)

}

Normalize <- function(boostRVC) {

  oldOrder <- diag(boostRVC$Matrix)
  Matrix <- reorderRVineMatrix(boostRVC$Matrix)

  out <- list(Matrix = Matrix,
              family = boostRVC$family,
              names = rev(boostRVC$names[oldOrder]),
              covariates = boostRVC$covariates)
  class(out) <- "boostRVineCopula"

  return(out)

}
