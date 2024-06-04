#' Calculate true log-likelihood of an Artificial Conditional R-Vine Copula Model
#'
#' @param U a matrix with d columns with values in [0,1].
#' @param response vector containing the linear predictor assumed among all bivariate conditional copulas.
#' @param boostRVM An \code{\link{boostRVineMatrix}} object including the R-vine structure in \code{Matrix},
#' the specified pair-copula families in \code{family} and the corresponding formulas in \code{formula}.
#' @param cores integer; the number of cores used for computations. Default: \code{cores = 1}.
#'
#' @return An \code{boostRVineCopula} object.
#'
#' @note The code of this function is mainly based on the R-package \code{\link{VineCopula}} by Nagler et al. but slightly adapted to the needs of this package.
#'
#' @export
# The following codes are mainly based on the R-package VineCopula by Nagler et al.
boostRVineTrueLL <- function(boostRVM,
                             U,
                             response,
                             cores = 1) {

  if (nrow(U) != length(response)) {
    stop("Rows of 'U' and length of 'response' need to be equal!")
  }

  d <- ncol(U)
  N <- nrow(U)

  ## reorder matrix to natural order
  family <- boostRVM$family
  M <- boostRVM$Matrix
  Mold <- M
  o <- diag(M)
  M <- reorderRVineMatrix(M)
  U_out <- U
  U <- U[, o[length(o):1]]

  ## create matrices required for selection of h-functions
  MaxMat <- createMaxMat(M)
  CondDistr <- neededCondDistr(M)

  ## create objects for results
  V <- list()
  V$direct <- array(NA, dim = c(d, N))
  V$indirect <- array(NA, dim = c(d, N))
  V$direct <- t(U[, d:1])
  pair_copulas <- lapply(1:d, function(k) lapply(1:d, function(j) NULL))

  for (k in d:2) {

    doEst <- function(i) {

      ## get pseudo-observaions
      m <- MaxMat[k, i]
      zr1 <- V$direct[i, ]

      zr2 <- if (m == M[k, i]) {
        V$direct[(d - m + 1), ]
      } else {
        V$indirect[(d - m + 1), ]
      }

      # parameter
      par1 <- tau2par(family[k, i], tanh(response))
      par2 <- 0
      # family
      fams <- getFams(as.numeric(family[k, i]))
      fam <- rep(fams[1], length(par1))
      if (length(fams) != 1) {
        fam[par1 < 0] <- fams[2]
      }

      ll <- log(BiCopPDF(u1 = zr2,
                         u2 = zr1,
                         family = fam,
                         par = par1,
                         par2 = par2,
                         check.pars = FALSE))



      direct <- indirect <- NULL
      if (CondDistr$direct[k - 1, i]) {

        direct <- suppressWarnings(BiCopHfunc1(u1 = zr2,
                                               u2 = zr1,
                                               family = fam,
                                               par = par1,
                                               par2 = par2,
                                               check.pars = FALSE))
      }


      if (CondDistr$indirect[k - 1, i]) {

        indirect <- suppressWarnings(BiCopHfunc2(u1 = zr2,
                                                 u2 = zr1,
                                                 family = fam,
                                                 par = par1,
                                                 par2 = par2,
                                                 check.pars = FALSE))
      }

      ## return results
      list(direct = direct,
           indirect = indirect,
           ll = ll)

    }
    ## run pair-copula selection for tree k
    res.k <- mclapply(seq_len(k - 1), doEst, mc.cores = cores)
    for (i in seq_len(k - 1)) {
      ## store info about selected pair-copula in list in same order as VineCopula does for matrices
      pair_copulas[[k]][[i]] <- res.k[[i]]$ll
      ## replace pseudo observations for estimation of next tree
      if (!is.null(res.k[[i]]$direct))
        V$direct[i, ] <- res.k[[i]]$direct
      if (!is.null(res.k[[i]]$indirect))
        V$indirect[i, ] <- res.k[[i]]$indirect
    } # end i = 1:(d-1)


    rm(res.k)
    gc()

  } # end k = d:2


  out <- pair_copulas

  return(out)

}


