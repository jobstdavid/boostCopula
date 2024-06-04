#' Sequential Gradient-Boosted Estimation of GLMs for a Conditional R-Vine Copula Model
#'
#' This function sequentially estimates the GLMs of conditional bivariate copulas of a d-dimensional conditional R-vine copula model
#' as specified by the corresponding \code{\link{boostRVineMatrix}} object via gradient-boosting.
#'
#' @param U a matrix with d columns with values in [0,1].
#' @param X matrix containing the covariates and having the same number of rows as \code{U}.
#' @param boostRVM An \code{\link{boostRVineMatrix}} object including the R-vine structure in \code{Matrix},
#' the specified pair-copula families in \code{family} and the corresponding formulas in \code{formula}.
#' @param control control parameters for the gradient-boosting estimation.
#' The \code{\link{boost_control}} function provides the default setting.
#' @param na.action a function which indicates what should happen if the data contains \code{NA}s. Default: \code{na.omit}.
#' @param light logical; if \code{FALSE} (default) the model frame is kept in the \code{boostBiCop} objects.
#' Otherwise, if \code{TRUE}, the model frame is omitted.
#' @param cores integer; the number of cores used for computations. Default: \code{cores = 1}.
#' @param ... unused.
#'
#' @return An \code{boostRVineCopula} object.
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
#' summary(object)
#'
#' @note The code of this function is mainly based on the R-package \code{\link{VineCopula}} by Nagler et al. but slightly adapted to the needs of this package.
#'
#' @importFrom stats AIC logLik
#' @importFrom VineCopula BiCopHfunc1 BiCopHfunc2
#' @export
# The following codes are mainly based on the R-package VineCopula by Nagler et al.
boostRVineSeqEst <- function(U,
                             X,
                             boostRVM,
                             control = boost_control(),
                             na.action = na.omit,
                             light = FALSE,
                             cores = 1,
                             ...) {

  if (nrow(U) != nrow(X)) {
    stop("Rows of U and X do not have equal length!")
  }

  # prepare na.action
  data <- cbind(U, X)
  data <- na.action(data)
  U <- data[, 1:ncol(U)]
  X <- data[, -c(1:ncol(U))]
  rm(data)

  d <- ncol(U)
  N <- nrow(U)
  ## set variable names
  if (is.null(colnames(U))) {
    colnames(U) <- paste("U", 1:ncol(U), sep = "")
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste("X", 1:ncol(X), sep = "")
  }

  ## reorder matrix to natural order
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

      cfit <- boostBiCopEst(formula = as.formula(boostRVM$formula[k, i]),
                            U = data.frame(U1 = zr2, U2 = zr1),
                            X = X,
                            family = boostRVM$family[k, i],
                            control = control,
                            na.action = na.pass,
                            light = light)

      ## transform data to pseudo-observations in next tree
      direct <- indirect <- NULL
      # parameter
      par1 <- predict(cfit, newdata = X, type = "parameter")$parameter
      par2 <- cfit$par2
      # family
      family <- getFams(as.numeric(cfit$family))
      fam <- rep(family[1], length(par1))
      if (length(family) != 1) {
        fam[par1 < 0] <- family[2]
      }

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
           cfit = cfit)

    }
    ## run pair-copula selection for tree k
    res.k <- mclapply(seq_len(k - 1), doEst, mc.cores = cores)
    for (i in seq_len(k - 1)) {
      ## store info about selected pair-copula in list in same order as VineCopula does for matrices
      pair_copulas[[k]][[i]] <- res.k[[i]]$cfit
      ## replace pseudo observations for estimation of next tree
      if (!is.null(res.k[[i]]$direct))
        V$direct[i, ] <- res.k[[i]]$direct
      if (!is.null(res.k[[i]]$indirect))
        V$indirect[i, ] <- res.k[[i]]$indirect
    } # end i = 1:(d-1)


    rm(res.k)
    gc()

  } # end k = d:2

  # calculate log-likelihood of model
  l <- unlist(lapply(1:d, function(k) Filter(Negate(is.null), pair_copulas[[k]])), recursive = F)
  loglik <- sum(sapply(l, logLik))
  df <- sum(sapply(1:(d*(d-1)/2), function(k) {attr(logLik(l[[k]]), "df")}))
  aic <- sum(sapply(l, AIC))

  out <- list(model_frame = data.frame(U_out, X),
              pair_copulas = pair_copulas,
              boostRVM = boostRVM,
              control = control,
              stats = list(n = N,
                           df = df,
                           loglik = loglik,
                           aic = aic))

  class(out) <- "boostRVineCopula"

  return(out)

}


###########################################################
# Helper functions for boostRVineSeqEst/boostRVineCopSelect
###########################################################
utils::globalVariables(c("Matrix"))
# command above helps to suppress annoying R-check notes for missing global variables

reorderRVineMatrix <- function(Matrix, oldOrder = NULL) {

  if (length(oldOrder) == 0) {
    oldOrder <- diag(Matrix)
  }
  O <- apply(t(1:nrow(Matrix)), 2, "==", Matrix)

  for (i in 1:nrow(Matrix)) {
    Matrix[O[, oldOrder[i]]] <- nrow(Matrix) - i + 1
  }
  return(Matrix)
}
createMaxMat <- function(Matrix) {

  if (dim(Matrix)[1] != dim(Matrix)[2])
    stop("Structure matrix has to be quadratic.")

  MaxMat <- reorderRVineMatrix(Matrix)

  n <- nrow(MaxMat)

  for (j in 1:(n - 1)) {
    for (i in (n - 1):j) {
      MaxMat[i, j] <- max(MaxMat[i:(i + 1), j])
    }
  }

  tMaxMat <- MaxMat
  tMaxMat[is.na(tMaxMat)] <- 0

  oldSort <- diag(Matrix)
  oldSort <- oldSort[n:1]

  for (i in 1:n) {
    MaxMat[tMaxMat == i] <- oldSort[i]
  }

  return(MaxMat)
}
neededCondDistr <- function(Vine) {

  if (dim(Vine)[1] != dim(Vine)[2])
    stop("Structure matrix has to be quadratic.")

  Vine <- reorderRVineMatrix(Vine)

  MaxMat <- createMaxMat(Vine)

  d <- nrow(Vine)

  M <- list()
  M$direct <- matrix(FALSE, d, d)
  M$indirect <- matrix(FALSE, d, d)

  M$direct[2:d, 1] <- TRUE

  for (i in 2:(d - 1)) {
    v <- d - i + 1

    bw <- as.matrix(MaxMat[i:d, 1:(i - 1)]) == v

    direct <- Vine[i:d, 1:(i - 1)] == v

    M$indirect[i:d, i] <- apply(as.matrix(bw & (!direct)), 1, any)

    M$direct[i:d, i] <- TRUE

    M$direct[i, i] <- any(as.matrix(bw)[1, ] & as.matrix(direct)[1, ])
  }

  return(M)
}

