#' Selection and Gradient-Boosted Estimation of GLMs for Conditional Bivariate Copulas
#'
#' Selection and gradient-boosted estimation of GLMs for conditional bivariate copulas.
#'
#' @param formula  an formula object, starting with the
#' \code{~} operator and followed by the covariates, which are separated by \code{+} operators.
#' @param U a matrix with two columns with values in [0,1].
#' @param X matrix containing the covariates and having the same number of rows as \code{U}.
#' @param familyset vector of bivariate copula families as described in \code{\link{boostBiCopEst}} to select from.
#' If \code{NA} (default), selection among all implemented copula families is performed.
#' @param selectioncrit character; either "\code{loglik}" or "\code{aic}" (default) for the copula family selection.
#' @param indeptest test for independence of both columns in \code{U}. Default: \code{NA}, i.e no independence check.
#' Otherwise, numeric value specifies significance level. The independence copula is chosen if the null hypothesis of independence cannot be rejected.
#' @param control control parameters for the gradient-boosting estimation.
#' The \code{\link{boost_control}} function provides the default setting.
#' @param na.action a function which indicates what should happen if the data contains \code{NA}s. Default: \code{na.omit}.
#' @param light logical; if \code{FALSE} (default) the model frame is kept in the \code{boostBiCop} object.
#' Otherwise, if \code{TRUE}, the model frame is omitted.
#' @param cores integer; the number of cores used for computations. Default: \code{cores = 1}.
#' @param ... unused.
#'
#' @return An \code{boostBiCop} object.
#'
#' @importFrom stats formula na.action na.omit optimize as.formula na.pass terms
#' @importFrom VineCopula BiCopIndTest
#' @importFrom utils tail
#' @importFrom Rfast colsums rowsums colmeans eachrow
#' @importFrom parallel mclapply
#'
#' @examples
#' # load simulated data
#' data(data_bicop)
#'
#' # fit object
#' object <- boostBiCopSelect(formula = ~.,
#'                            U = data_bicop[, 1:2],
#'                            X = data_bicop[, -c(1:2)],
#'                            familyset = NA,
#'                            control = boost_control(deselection = "attributable"))
#' summary(object)
#' plot(object)
#'
#' @export
boostBiCopSelect <- function(formula,
                             U,
                             X,
                             familyset = NA,
                             selectioncrit = "aic",
                             indeptest = NA,
                             control = boost_control(),
                             na.action = na.omit,
                             light = FALSE,
                             cores = 1,
                             ...) {
  # prepare data
  data <- data_prepare(formula = formula,
                       U = U,
                       X = X,
                       na.action = na.action)
  # center data
  Z <- data_center(data = data,
                   center = control$center)

  # check for Independence copula
  p.value <- BiCopIndTest(Z$u[, 1], Z$u[, 2])$p.value
  # create Independence copula
  if (!is.na(indeptest) & (p.value >= indeptest)) {

    coefs <- rep(0, ncol(Z$x))
    names(coefs) <- colnames(Z$x)
    object <- list(family = "0",
                   familyname = "Independence",
                   formula = formula,
                   control = control,
                   location = Z$location,
                   coefficients = coefs,
                   par2 = 0,
                   iterations = 0,
                   paths = list(coef_path = matrix(coefs, nrow = 1),
                                loglik_path = 0,
                                cov_path = NA),
                   stats = list(n = nrow(Z$x),
                                df = 0,
                                loglik = 0,
                                aic = 0)
                   )
    if (!light) {
      object <- append(object,
                       list(model_frame = data),
                           after = 3)
    }
    class(object) <- "boostBiCop"

  } else { # fit copula families 1, 301:304, 401:404

    if (is.na(familyset[1])) {
      familyset <- c(1, 301:304, 401:404)
    }
    d <- length(familyset)
    familyset <- lapply(1:d, function(j) fam_trafo(family = familyset[[j]]))

    m <- mclapply(1:d, function(j) boostBiCop(u = Z$u, x = Z$x, family = familyset[[j]], control = control), mc.cores = cores)

    if (control$mstop == "aic") {

      # find optimal number of stopping iterations by aic for each family
      s_m <- lapply(1:d, function(j) aic(m[[j]], control))
      # deselection procedure
      if (control$deselection %in% c("attributable", "cumulative")) {
        m <- lapply(1:d, function(j) subset_model(m[[j]], s_m[[j]]$iterations))
        m <- mclapply(1:d, function(j) desel_intern(m[[j]],
                                                    Z$u,
                                                    Z$x,
                                                    familyset[[j]],
                                                    control),
                      mc.cores = cores)
        # get scores for copula family selection
        s_m <- lapply(1:d, function(j) {
          iterations <- length(m[[j]]$paths$cov_path)

          list(iterations = iterations,
               loglik = m[[j]]$paths$loglik_path[iterations],
               aic = -2*m[[j]]$paths$loglik_path[iterations]+2*sum(m[[j]]$paths$coef_path[iterations, ] != 0))

        })
      }
      # get scores for copula family selection
      sel.score <- sapply(1:d, function(j) selcrit.score(s_m[[j]], selectioncrit))
      # get optimal copula family
      index <- which.max(sel.score)
      # get optimal number of stopping iterations for selected copula family
      iterations <- s_m[[index]]$iterations
      # select copula family using optimal number of stopping iterations
      m <- subset_model(m[[index]], iterations)

    } else if (control$mstop == "cv") {

      # find optimal number of stopping iterations by cv for each family
      s_m <- mclapply(1:d, function(j) cv(m[[j]], Z$u, Z$x, familyset[[j]], control), mc.cores = cores)
      # deselection procedure
      if (control$deselection %in% c("attributable", "cumulative")) {
        m <- lapply(1:d, function(j) subset_model(m[[j]], s_m[[j]]$iterations))
        m <- mclapply(1:d, function(j) desel_intern(m[[j]],
                                                    Z$u,
                                                    Z$x,
                                                    familyset[[j]],
                                                    control),
                      mc.cores = cores)
        # get scores for copula family selection
        s_m <- lapply(1:d, function(j) {
          iterations <- length(m[[j]]$paths$cov_path)

          list(iterations = iterations,
               loglik = m[[j]]$paths$loglik_path[iterations],
               aic = -2*m[[j]]$paths$loglik_path[iterations]+2*sum(m[[j]]$paths$coef_path[iterations, ] != 0))

        })
      }
      # get scores for copula family selection
      sel.score <- sapply(1:d, function(j) selcrit.score(s_m[[j]], selectioncrit))
      # get optimal copula family
      index <- which.max(sel.score)
      # get optimal number of stopping iterations for selected copula family
      iterations <- s_m[[index]]$iterations
      # select copula family using optimal number of stopping iterations
      m <- subset_model(m[[index]], iterations)

    } else if (control$mstop == "max") {

      # maximal number of iterations
      iterations <- control$maxit
      # deselection procedure
      if (control$deselection %in% c("attributable", "cumulative")) {
        m <- lapply(1:d, function(j) subset_model(m[[j]], iterations))
        m <- mclapply(1:d, function(j) desel_intern(m[[j]],
                                                    Z$u,
                                                    Z$x,
                                                    familyset[[j]],
                                                    control),
                      mc.cores = cores)
      }
      # get scores for copula family selection
      s_m <- lapply(1:d, function(j) {list(loglik = m[[j]]$paths$loglik_path[iterations],
                                           aic = -2*m[[j]]$paths$loglik_path[iterations]+2*sum(m[[j]]$paths$coef_path[iterations, ] != 0))})
      sel.score <- sapply(1:d, function(j) selcrit.score(s_m[[j]], selectioncrit))
      # get optimal copula family
      index <- which.max(sel.score)
      # select copula family using optimal number of stopping iterations
      m <- subset_model(m[[index]], iterations)

    } else {
      stop("This stopping criterion is not available!")
    }

    # coefficients are the same for model with centered and non-centered covariates
    # except of the intercept which needs to be adapted for the centered-covariates case to obtain the correct non-centered model
    if (control$center) {
      if(attr(terms(formula, allowDotAsName = T), "intercept") != 0) { # contains intercept
        m$coefficients[1] <- m$coefficients[1] - sum(m$coefficients[-1]*Z$location[names(m$coefficients[-1])])
        m$paths$coef_path[, 1] <- m$paths$coef_path[, 1] - matvecmult_eigen(matrix(m$paths$coef_path[, -1], ncol = length(names(m$coefficients[-1]))), Z$location[names(m$coefficients[-1])], control$cores.control)
      }
    }

    # add initialization step 0
    m$paths$coef_path <- rbind(0, m$paths$coef_path)
    m$paths$cov_path <- c(NA, m$paths$cov_path)
    m$paths$loglik_path <- c(-sum(familyset[[index]]$loss(u = Z$u, f = rep(0, nrow(Z$u)), par2 = familyset[[index]]$par2)),
                             m$paths$loglik_path)

    object <- finalize_boostBiCop(formula = formula,
                                  m = m,
                                  data = data,
                                  control = control,
                                  location = Z$location,
                                  light = light)

  }

  return(object)


}



# function for data preparation
data_prepare <- function(formula, U, X, na.action)  {

  # check U
  # if (class(U) != "data.frame") {
  #   stop("U must be a data frame!")
  # }

  if (ncol(U) != 2) {
    stop("U must have two columns!")
  }

  if(is.null(colnames(U))) {
    colnames(U) <- paste0("U", 1:2)
  }

  # check X
  # if (class(X) != "data.frame") {
  #   stop("X must be a data frame!")
  # }

  if (nrow(U) != nrow(X)) {
    stop("Rows of U and X do not have equal length!")
  }

  if(is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:ncol(X))
  }

  # check model matrix (faster than model.matrix)
  if (formula != ~.) {
    covs <- all.vars(formula)
    if (any(!covs %in% colnames(X))) {
      stop("Covariates in X do not match with formula!")
    }
    X <- as.matrix(X[, covs])
    colnames(X) <- covs
    if (attr(terms(formula), "intercept") != 0) {
      X <- cbind(1, X)
      colnames(X)[1] <- "(Intercept)"
    }
  } else {
    X <- cbind(1, X)
    colnames(X)[1] <- "(Intercept)"
  }

  # combine, apply na.action and return data
  return(na.action(cbind(U, X)))

}
# function for data centering
data_center <- function(data, center) {

  # get pseudo-observations
  U <- data[, c(1:2)]
  X_names <- colnames(data)[-c(1,2)]

  # (not) center covariates
  if (center) {
    intercept <- which(colnames(data) == "(Intercept)")
    X <- as.matrix(data[, -c(1, 2, intercept)])

    location <- colmeans(X)
    names(location) <- colnames(data)[-c(1, 2, intercept)]
    X <- eachrow(X, location, oper = "-")

    if (length(intercept) > 0) {
      X <- cbind(1, X)
    }
  } else {
    X <- as.matrix(data[, -c(1,2)])
    location <- NULL
  }

  colnames(X) <- X_names
  out <- list(u = U,
              x = X,
              location = location)
  return(out)
}
# function to obtain the copula family
fam_trafo <- function(family) {

  if (family == 1) {
    BiCop_1_tau()
  } else if (family == 301) {
    BiCop_301_tau()
  } else if (family == 302) {
    BiCop_302_tau()
  } else if (family == 303) {
    BiCop_303_tau()
  } else if (family == 304) {
    BiCop_304_tau()
  } else if (family == 401) {
    BiCop_401_tau()
  } else if (family == 402) {
    BiCop_402_tau()
  } else if (family == 403) {
    BiCop_403_tau()
  } else if (family == 404) {
    BiCop_404_tau()
  } else {
    stop("This family is not implemented!")
  }

}

# function for gradient-boosting
boostBiCop <- function(u, x, family, control) {

  # initialize copula family
  ngr <- family$ngradient
  loss <- family$loss
  par2 <- family$par2

  # initialize settings
  maxit <- control$maxit
  nu <- control$nu
  n_cores <- control$cores.control

  # constant for boosting: beta_j = <ngrad,X_j>^2/C^2
  C <- sqrt(colsums(x^2))

  # initialize parameters
  beta <- rep(0, ncol(x))
  f <- rep(0, nrow(x))
  loglik_path <- coef_path <- c()
  cov_path <- c()
  fill <- FALSE

  # start boosting
  for (k in 1:maxit) {

    ngrad <- ngr(u = u, f = f, par2 = par2)
    boost_cor <- matTvecmult_eigen(x, ngrad, n_cores)/C
    # stop iteration, if infinite values arise
    if (any(c(!is.finite(ngrad), !is.finite(boost_cor), !is.finite(loglik_path)))) {
      fill <- TRUE
      break
    }
    index <- which.max(abs(boost_cor))
    bl <- nu*boost_cor[index]/C[index]
    beta[index] <- beta[index] + bl
    f <- f + bl*x[, index]

    # track boosting
    cov_path <- c(cov_path, index)
    coef_path <- rbind(coef_path, beta)
    loglik_path <- c(loglik_path, -sum(loss(u = u, f = f, par2 = par2)))

  }

  if (fill) {
    cov_path <- c(cov_path, rep(NA, maxit-length(cov_path)))
    coef_path <- rbind(coef_path, matrix(rep(beta, maxit-length(loglik_path)), ncol = length(beta), byrow = T))
    loglik_path <- c(loglik_path, rep(tail(loglik_path, 1), maxit-length(loglik_path)))
  }

  colnames(coef_path) <- colnames(x)
  rownames(coef_path) <- NULL
  n <- nrow(u)

  out <- list(family = family$family,
              familyname = family$familyname,
              paths = list(coef_path = coef_path,
                           loglik_path = loglik_path,
                           cov_path = cov_path),
              par2 = par2,
              n = n)

  return(out)

}

# function for optimal number of iterations based on AIC; returns iterations, AIC, loglik
aic <- function(m, control) {

  # dof
  df <- rowsums(m$paths$coef_path != 0)
  # AIC
  aic <- -2*m$paths$loglik_path+2*df
  # find stopping-iteration
  iterations <- which.min(aic)
  loglik <- m$paths$loglik_path[iterations]
  aic <- -2*loglik + 2*df[iterations]

  out <- list(iterations = iterations,
              loglik = loglik, # loglik based on iterations
              aic = aic) # aic based on iterations

  return(out)

}
# function for calculation the score of the copula family selection criterion
selcrit.score <- function(opt, selectioncrit) {

  if (selectioncrit == "loglik") {
    opt$loglik
  } else if (selectioncrit == "aic") {
    -opt$aic # minus in front of opt$aic, to be able to apply which.max afterwards
  } else {
    stop("This 'selectioncrit' is not implemented!")
  }

}
# function for model subsetting
subset_model <- function(m, iterations) {

  # obtain model based on stopping iterations
  coefficients <- m$paths$coef_path[iterations, ]
  cov_path <- m$paths$cov_path[1:iterations]
  loglik_path <- m$paths$loglik_path[1:iterations]
  coef_path <- matrix(m$paths$coef_path[1:iterations, ],
                      nrow = iterations,
                      dimnames = list(c(), colnames(m$paths$coef_path)))
  colnames(coef_path) <- colnames(m$paths$coef_path)
  loglik <- tail(loglik_path, 1)
  df <- sum(coefficients != 0)
  aic <- -2*loglik+2*df
  family <- m$family
  familyname <- m$familyname
  par2 <- m$par2
  n <- m$n

  out <- list(family = family,
              familyname = familyname,
              coefficients = coefficients,
              par2 = par2,
              iterations = iterations,
              paths = list(coef_path = coef_path,
                           loglik_path = loglik_path,
                           cov_path = cov_path),
              stats = list(n = n,
                           df = df,
                           loglik = loglik,
                           aic = aic))


  return(out)

}
# function for optimal number of iterations based on CV
cv <- function(m, u, x, family, control) {

  # initial settings
  foldid <- control$foldid
  nfolds <- control$nfolds
  # in case gradient-boosting stopped earlier due to negative values
  # coefficients stay constant after this stopping iteration => no improvement
  # therefore set the optimal number of iterations to the minimum of
  # the earlier stopping iteration
  if (any(is.na(m$paths$cov_path))) {
    control$maxit <- sum(!is.na(m$paths$cov_path))
  }
  maxit <- control$maxit

  if (maxit > 1) {

    mc.cores <- control$cores.control
    d <- ncol(x)
    par2 <- family$par2

    if(is.null(foldid)) {
      foldid <- sample(1:nfolds, size = m$n, replace = TRUE)
    } else {
      nfolds <- length(unique(foldid))
    }

    # start cv
    loglik_test <- mclapply(1:nfolds, function(i) {

      # indices (rows) training and test data set
      train <- foldid != unique(foldid)[i]
      test <- !train

      # cv-model
      m_cv <- boostBiCop(u[train, ], as.matrix(x[train, ]), family, control)
      # copula parameter vector for test data
      f_test <- as.vector(matmatTmult_eigen(as.matrix(x[test, ]), m_cv$paths$coef_path, 1))
      # pseudo-copula observations adapted to f_test
      u_test <- matrix(rep(t(u[test, ]), maxit), ncol = 2, byrow = TRUE)
      # negative log-likelihood
      loglik <- -family$loss(u_test, f_test, par2)
      # negative log-likelihood matrix for ith-fold
      matrix(loglik, ncol = maxit, byrow = FALSE)

    }, mc.cores = mc.cores)
    # end cv

    # find stopping-iteration
    iterations <- which.max(colsums(do.call(rbind, loglik_test)))
    if (length(iterations) == 0) {
      iterations <- maxit
    }
  } else {
    iterations <- maxit
  }
  # get logLik based on stopping-iterations
  loglik <- m$paths$loglik_path[iterations]
  aic <- -2*loglik + 2*sum(m$paths$coef_path[iterations, ] != 0)

  out <- list(iterations = iterations,
              loglik = loglik, # loglik based on iterations
              aic = aic) # aic based on iterations

  return(out)

}

# finalize boostBiCop object
finalize_boostBiCop <- function(formula, m, data, control, location, light) {

  covs <- names(m$coefficients)
  if (control$deselection %in% c("attributable", "cumulative")) {
    if (length(covs) == 1 & "(Intercept)" %in% covs) {
      formula <- ~ 1
    } else {
      formula <- as.formula(paste0("~ ", paste0(covs[covs != "(Intercept)"], collapse = "+")))
    }
  }


  m <- append(m,
              list(formula = formula,
                   control = control,
                   location = location[covs[covs != "(Intercept)"]]),
              after = 2)

  if (!light) {
    vars <- colnames(data)
    if (control$deselection %in% c("attributable", "cumulative")) {
      vars <- c(vars[1:2], names(m$coefficients))
    }
    m <- append(m,
                list(model_frame = data[, vars]),
                after = 3)
  }

  class(m) <- "boostBiCop"
  return(m)

}

# deselection function via "attributable" or "cumulative" risk
desel_intern <- function(m, u, x, family, control) {

  iterations <- m$iterations
  if (iterations > 1) {

    covs <- names(m$coefficients)
    loglik_path <- c(-sum(family$loss(u = u, f = rep(0, m$stats$n), par2 = family$par2)), m$paths$loglik_path)
    total_loglik <- loglik_path[iterations+1] - loglik_path[1]
    loglik_diff <- diff(loglik_path)

    select <- m$paths$cov_path - 1

    # remove risk reduction for intercept
    loglik_diff <- loglik_diff[select != 0]
    # remove intercept positions
    select <- select[select != 0]

    # sort selected variables
    covs <- sort(unique(select))
    # calculate attributable risk reduction for each variable
    loglik.covs <- sapply(1:length(covs), function(j) sum(loglik_diff[which(covs[j] == select)]) )

    # remove intercept from covariates
    pars <- names(m$coefficients[m$coefficients != 0])
    if ("(Intercept)" %in% pars) {
      pars <- pars[-which(pars == "(Intercept)")]
    }

    if (length(covs) != 0) {

      loglik.order <- data.frame(covs, pars, loglik.covs)
      loglik.order <- loglik.order[order(loglik.order$loglik.covs, decreasing = T), ]
      loglik.order$cumloglik <- cumsum(loglik.order$loglik.covs)
      colnames(loglik.order) <- c("cov", "covname", "loglik", "cumloglik")

      gamma <- control$gamma
      gamma_loglik <- total_loglik * gamma

      method <- control$deselection
      if (method == "attributable"){
        delete <- which(loglik.order$loglik < gamma_loglik)
        if (length(delete) > 0) {
          sel <- loglik.order[-delete, ]
        } else {
          sel <- loglik.order
        }
      } else if (method == "cumulative") {
        delete <- which(loglik.order$cumloglik < gamma_loglik)
        if (length(delete) > 0) {
          sel <- loglik.order[-delete, ]
        } else {
          sel <- loglik.order
        }
      }

      if (nrow(sel) == 0) {
        x_sel <- "(Intercept)"
      } else{
        sel.cov <- as.character(sel$covname)
        covs <- names(m$coefficients)
        x_sel <- covs[covs %in% c("(Intercept)", sel.cov)]
      }

      x <- as.matrix(x[, x_sel])
      colnames(x) <- x_sel
      control$maxit <- m$iterations

      m <- boostBiCop(u = u, x = x, family = family, control = control)

    }

  }

  return(m)

}


