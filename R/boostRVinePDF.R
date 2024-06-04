#' Density of a Conditional Regular Vine Copula
#'
#' Probability density function (PDF) of an \code{boostRVineCopula} object.
#'
#' @param object an \code{boostRVineCopula} object.
#' @param U a matrix with two columns with values in [0,1].
#' @param X matrix containing the covariates and having the same number of rows as \code{U}.
#' @param cores integer; the number of cores used for computations. Default: \code{cores = 1}.
#' @param ... unused.
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
#' boostRVinePDF(object,
#'               U = data_vinecop[, 1:5],
#'               X = data_vinecop[, -c(1:5)])
#'
#' @return A numeric vector containing the density values of the considered \code{boostRVineCopula} object, evaluated at \code{U} and \code{X}.
#'
#' @note The code of this function is mainly based on the R-package \code{\link{VineCopula}} by Nagler et al. but slightly adapted to the needs of this package.
#'
#' @importFrom stats predict
#' @export

boostRVinePDF <- function(object, U, X, cores = 1, ...) {

  if(missing(U) | missing(X)) {
    if (!"model_frame" %in% names(object)) {
      stop("Missing 'model_frame' in object due to 'light = TRUE'! Please provide arguments for 'U' and 'X'!")
    }
    U <- object$model_frame[, 1:ncol(object$boostRVM$Matrix)]
    X <- object$model_frame[, -c(1:ncol(object$boostRVM$Matrix))]
  }

  if (ncol(U) != ncol(object$boostRVM$Matrix)) {
    stop(paste0("U must have ", ncol(object$boostRVM$Matrix), " two columns!"))
  }

  if (nrow(U) != nrow(X)) {
    stop("Rows of U and X do not have equal length!")
  }

  # prepare na.action
  data <- cbind(U, X)
  data <- na.omit(data)
  U <- data[, 1:ncol(U)]
  X <- data[, -c(1:ncol(U))]
  rm(data)

  d <- ncol(U)
  N <- nrow(U)

  ## reorder matrix to natural order
  M <- object$boostRVM$Matrix
  Mold <- M
  o <- diag(M)
  M <- reorderRVineMatrix(M)
  U <- U[, o[length(o):1]]

  ## create matrices required for selection of h-functions
  MaxMat <- createMaxMat(M)
  CondDistr <- neededCondDistr(M)

  ## create objects for results
  V <- list()
  V$direct <- array(NA, dim = c(d, N))
  V$indirect <- array(NA, dim = c(d, N))
  V$direct <- t(U[, d:1])
  pc <- object$pair_copulas
  density_mat <- c()

  for (k in d:2) {

    pdf_tree <- function(i) {

      ## get pseudo-observations
      m <- MaxMat[k, i]
      zr1 <- V$direct[i, ]

      zr2 <- if (m == M[k, i]) {
        V$direct[(d - m + 1), ]
      } else {
        V$indirect[(d - m + 1), ]
      }

      cfit <- pc[[k]][[i]]

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

      density <- BiCopPDF(u1 = zr2,
                          u2 = zr1,
                          family = fam,
                          par = par1,
                          par2 = par2,
                          check.pars = FALSE)

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
           density = density)

    }
    ## run pair-copula selection for tree k
    res.k <- mclapply(seq_len(k - 1), pdf_tree, mc.cores = cores)

    for (i in seq_len(k - 1)) {
      ## store info about selected pair-copula in list in same order as VineCopula does for matrices
      density_mat <- c(density_mat, res.k[[i]]$density)
      ## replace pseudo observations for estimation of next tree
      if (!is.null(res.k[[i]]$direct))
        V$direct[i, ] <- res.k[[i]]$direct
      if (!is.null(res.k[[i]]$indirect))
        V$indirect[i, ] <- res.k[[i]]$indirect
    } # end i = 1:(d-1)


    rm(res.k)
    gc()

  } # end k = d:2

  # calculate pdf values of model
  density_mat <- matrix(density_mat, nrow = N, ncol = d*(d-1)/2, byrow = FALSE)
  out <- apply(density_mat, 1, prod)

  return(out)

}
