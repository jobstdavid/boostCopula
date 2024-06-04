#' Sequential Selection and Gradient-Boosted Estimation of GLMs for a Conditional R-Vine Copula Model
#'
#' This function sequentially estimates GLMs for conditional bivariate copulas of a d-dimensional conditional R-vine copula model
#' as specified by the corresponding \code{\link{boostRVineMatrix}} object via gradient-boosting.
#' Additionally, it selects the best fitting copula family for each conditional bivariate copula.
#'
#' @param U a matrix with d columns with values in [0,1].
#' @param X matrix containing the covariates and having the same number of rows as \code{U}.
#' @param boostRVM An \code{\link{boostRVineMatrix}} object.
#' Note that the \code{family} matrix in the \code{boostRVineMatrix} object will be ignored to select from the \code{familyset}.
#' @param familyset vector of bivariate copula families as described in \code{\link{boostBiCopEst}} to select from. If \code{NA} (default),
#' selection among all implemented copula families is performed. Note, that a defined \code{family} matrix in \code{boostRVM} will be ignored
#' and therefore does not need to be specified.
#' @param selectioncrit character; either "\code{loglik}" or "\code{aic}" (default) for the copula family selection.
#' @param indeptest test for independence for each edge in the regular vine. Default: \code{NA}, i.e no independence check.
#' Otherwise, numeric value specifies significance level. The independence copula is chosen if the null hypothesis of independence cannot be rejected.
#' @param trunclevel integer; level of truncation.
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
#' Matrix <- matrix(Matrix, 5, 5, byrow = TRUE)
#' Family <- matrix(0, nrow = 5, ncol = 5)
#' Formula <- matrix("~ .", 5, 5, byrow = TRUE)
#' boostRVM <- boostRVineMatrix(Matrix = Matrix,
#'                              family = Family,
#'                              formula = Formula)
#'
#' # fit object
#' object <- boostRVineCopSelect(U = data_vinecop[, 1:5],
#'                               X = data_vinecop[, -c(1:5)],
#'                               boostRVM = boostRVM,
#'                               familyset = c(1, 301, 304, 401, 404),
#'                               selectioncrit = "aic",
#'                               control = boost_control(deselection = "attributable"),
#'                               cores = 1)
#' summary(object)
#'
#'
#' @note The code of this function is mainly based on the R-package \code{\link{VineCopula}} by Nagler et al. but slightly adapted to the needs of this package.
#'
#' @export
# The following codes are mainly based on the R-package VineCopula by Nagler et al.
boostRVineCopSelect <- function(U,
                                X,
                                boostRVM,
                                familyset = NA,
                                selectioncrit = "aic",
                                indeptest = NA,
                                trunclevel = NA,
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

  # ignore 'family' in boostRVM, as this will be selected out of familyset
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

    if (!is.na(trunclevel) & trunclevel <= (d-k)) {
      indeptest <- -Inf
    }

    doEst <- function(i) {

      ## get pseudo-observations
      m <- MaxMat[k, i]
      zr1 <- V$direct[i, ]

      zr2 <- if (m == M[k, i]) {
        V$direct[(d - m + 1), ]
      } else {
        V$indirect[(d - m + 1), ]
      }

      # prepare data with covariates
      cfit <- boostBiCopSelect(as.formula(boostRVM$formula[k, i]),
                               U = cbind(zr2, zr1),
                               X = X,
                               familyset = familyset,
                               selectioncrit = selectioncrit,
                               indeptest = indeptest,
                               control = control,
                               na.action = na.pass,
                               light = light,
                               cores = cores)

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
    res.k <- lapply(seq_len(k - 1), function(i) {doEst(i)})
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




