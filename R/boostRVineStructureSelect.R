#' Sequential Gradient-Boosted Estimation of GLMs for a Conditional R-Vine Copula Model
#'
#' This function sequentially estimates the GLMs of conditional bivariate copulas as well as
#' the structure of a d-dimensional conditional R-vine copula model via gradient-boosting.
#' Additionally, it selects the best fitting copula family for each conditional bivariate copula.
#'
#' @param U a matrix with d columns with values in [0,1].
#' @param X matrix containing the covariates and having the same number of rows as \code{U}.
#' @param boostRVM An \code{\link{boostRVineMatrix}} object.
#' Note that the \code{family} in the \code{boostRVineMatrix} object will be ignored to select from the \code{familyset}.
#' Additionally, the regular vine structure \code{Matrix} will be ignored.
#' @param familyset vector of bivariate copula families as described in \code{\link{boostBiCopEst}} to select from. If \code{NA} (default),
#' selection among all implemented copula families is performed.
#' @param selectioncrit character; either "\code{loglik}" or "\code{aic}" (default) for the copula family selection.
#' @param indeptest test for independence for each edge in the regular vine. Default: \code{NA}, i.e no independence check.
#' Otherwise, numeric value specifies significance level. The independence copula is chosen if the null hypothesis of independence cannot be rejected.
#' @param trunclevel integer; level of truncation.
#' @param treecrit character; edge weight for Dissman's structure selection algorithm.
#' Either "\code{tau}" (default) for absolute value of empirical Kendall's tau or
#' "\code{rho}" for absolute value of empirical Spearmans's rho.
#' @param vine_type integer; type of the vine model to be specified:\cr
#' \code{0} for R-vine (default)\cr
#' \code{1} for C-vine.
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
#' Matrix <- matrix(0, 5, 5, byrow = TRUE)
#' Family <- matrix(0, nrow = 5, ncol = 5)
#' Formula <- matrix("~ .", 5, 5, byrow = TRUE)
#' boostRVM <- boostRVineMatrix(Matrix = Matrix,
#'                              family = Family,
#'                              formula = Formula)
#'
#' # fit object
#' object <- boostRVineStructureSelect(U = data_vinecop[, 1:5],
#'                                     X = data_vinecop[, -c(1:5)],
#'                                     boostRVM = boostRVM,
#'                                     familyset = c(1, 301, 304, 401, 404),
#'                                     selectioncrit = "aic",
#'                                     treecrit = "tau",
#'                                     vine_type = 0,
#'                                     control = boost_control(deselection = "attributable"),
#'                                     cores = 1)
#' summary(object)
#'
#'
#' @note The code of this function is mainly based on the R-package \code{\link{VineCopula}} by Nagler et al. but slightly adapted to the needs of this package.
#'
#' @importFrom stats cor
#' @importFrom utils combn
#' @export
# The following codes are mainly based on the R-package VineCopula by Nagler et al.
boostRVineStructureSelect <- function(U,
                                      X,
                                      boostRVM,
                                      familyset = NA,
                                      selectioncrit = "aic",
                                      indeptest = NA,
                                      trunclevel = NA,
                                      treecrit = "tau",
                                      vine_type = 0,
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
  varnames <- colnames(U)
  covnames <- colnames(X)


  ## sanity checks
  if (d <= 2)
    stop("Dimension has to be at least 3! Use function 'boostBiCopSelect' instead!")

  ## set defaults
  if (vine_type == 0) {
    vine_type <- "RVine"
  } else {
    vine_type <- "CVine"
  }
  if (is.na(trunclevel)) {
    trunclevel <- d
  }

  # save settings for final output
  settings <- list(boostRVM = boostRVM,
                   familyset = familyset,
                   selectioncrit = selectioncrit,
                   indeptest = indeptest,
                   trunclevel = trunclevel,
                   treecrit = treecrit,
                   vine_type = vine_type,
                   control = control,
                   na.action = na.action,
                   light = light)

  treecrit <- set_treecrit(treecrit)
  ## initialize object for results
  RVine <- list(Tree = NULL, Graph = NULL)

  ## estimation in first tree ----------------------------
  # find optimal tree
  g <- initializeFirstGraph(U, treecrit)
  MST <- findMaxTree(g, mode = vine_type)

  # estimate pair-copulas
  VineTree <- fit.FirstTreeCopulas(MST,
                                   U,
                                   X,
                                   boostRVM$formula[d, ],
                                   familyset,
                                   selectioncrit,
                                   indeptest,
                                   trunclevel,
                                   control,
                                   light,
                                   cores)


  RVine$Tree[[1]] <- VineTree
  RVine$Graph[[1]] <- g
  oldVineGraph <- VineTree

  ## estimation in first tree ----------------------------
  # find optimal tree
  g <- initializeFirstGraph(U, treecrit)
  MST <- findMaxTree(g, mode = vine_type)

  # estimate pair-copulas
  VineTree <- fit.FirstTreeCopulas(MST,
                                   U,
                                   X,
                                   boostRVM$formula[d, ],
                                   familyset,
                                   selectioncrit,
                                   indeptest,
                                   trunclevel,
                                   control,
                                   light,
                                   cores)


  RVine$Tree[[1]] <- VineTree
  RVine$Graph[[1]] <- g
  oldVineGraph <- VineTree

  ## estimation in higher trees --------------------------
  for (tree in 2:(d - 1)) {

    ## old pseudo-observations are uneccessary in RVine object
    RVine$Tree[[tree - 1]]$E$Copula.CondData.1 <- NULL
    RVine$Tree[[tree - 1]]$E$Copula.CondData.2 <- NULL

    if (!is.na(trunclevel) & trunclevel == (tree - 1)) {
      indeptest <- -Inf
    }

    # find optimal tree
    g <- buildNextGraph(VineTree, treecrit, truncated = trunclevel < tree)
    MST <- findMaxTree(g, mode = vine_type, truncated = trunclevel < tree)
    # estimate pair-copulas
    VineTree <- fit.TreeCopulas(X,
                                MST,
                                VineTree,
                                boostRVM$formula[d+1-tree, ],
                                familyset,
                                selectioncrit,
                                indeptest,
                                trunclevel,
                                control,
                                light,
                                cores)
    # store results
    RVine$Tree[[tree]] <- VineTree
    RVine$Graph[[tree]] <- g

  }

  gc()


  out <- as.boostRVineCopula(RVine, U, X, settings)

  return(out)

}


###########################################################
# Helper functions for boostRVineStructureSelect
###########################################################

# initializing functions
set_treecrit <- function(treecrit) {
  ## check if function is appropriate or type is implemented
  if (all(treecrit == "tau")) {

    treecrit <- function(u1, u2) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 10) {
        tau <- 0
      } else {
        complete.freq <- mean(!is.na(u1 + u2))
        tau <- abs(cor(u1[complete.i], u2[complete.i], method = "kendall"))
        tau <- tau * sqrt(complete.freq)
      }
      tau
    }

  } else if (all(treecrit == "rho")) {

    treecrit <- function(u1, u2) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 10) {
        rho <- 0
      } else {
        complete.freq <- mean(!is.na(u1 + u2))
        rho <- abs(cor(u1, u2, method = "spearman", use = "complete.obs"))
        rho <- rho * sqrt(complete.freq)
      }
      rho
    }

  }

  ## return treecrit function
  treecrit

}
initializeFirstGraph <- function(data, treecrit) {
  ## calculate edge weight for each possible pair
  all.pairs <- combn(1:ncol(data), 2)
  edge.ws <- apply(all.pairs, 2,
                   function(ind)
                     treecrit(data[, ind[1]], data[, ind[2]]))
  # number of pairwise complete observations / all observations
  rel.nobs <- apply(all.pairs, 2,
                    function(ind)
                      mean(!is.na(data[, ind[1]] + data[, ind[2]])))
  edge.ws <- edge.ws


  ## store in symmetric matrix with appropriate names
  W <- diag(ncol(data))
  W[lower.tri(W)] <- edge.ws
  W <- t(W)
  colnames(W) <- rownames(W) <- colnames(data)

  ## return full graph with edge weights
  graphFromWeightMatrix(W)
}
findMaxTree <- function(g, mode = "RVine", truncated = FALSE) {

  if (truncated == FALSE) {
    ## construct adjency matrix
    A <- adjacencyMatrix(g)
    d <- ncol(A)

    if (mode == "RVine") {
      ## initialize
      tree <- NULL
      edges <- matrix(NA, d - 1, 2)
      w <- numeric(d - 1)
      i <- 1

      ## construct minimum spanning tree
      for (k in 1:(d - 1)) {
        # add selected edge to tree
        tree <- c(tree, i)

        # find edge with minimal weight
        m <- apply(as.matrix(A[, tree]), 2, min)
        a <- apply(as.matrix(A[, tree]), 2,
                   function(x) order(rank(x)))[1, ]
        b <- order(rank(m))[1]
        j <- tree[b]
        i <- a[b]

        # store edge and weight
        edges[k, ] <- c(j, i)
        w[k] <- A[i, j]

        ## adjust adjecency matrix to prevent loops
        for (t in tree)
          A[i, t] <- A[t, i] <- Inf
      }

      ## reorder edges for backwads compatibility with igraph output
      edges <- t(apply(edges, 1, function(x) sort(x)))
      edges <- edges[order(edges[, 2], edges[, 1]), ]

      ## delete unused edges from graph
      E <- g$E$nums
      in.tree <- apply(matrix(edges, ncol = 2), 1,
                       function(x) which((x[1] == E[, 1]) & (x[2] == E[, 2])))
      MST <- g
      g$E$todel <- rep(TRUE, nrow(E))
      if (any(g$E$todel)) {
        g$E$todel[in.tree] <- FALSE
        MST <- deleteEdges(g)
      }
    } else if (mode  == "CVine") {
      ## set root as vertex with minimal sum of weights
      A <- adjacencyMatrix(g)
      diag(A) <- 0
      sumtaus <- rowSums(A)
      root <- which.min(sumtaus)

      ## delete unused edges
      g$E$todel <- !((g$E$nums[, 2] == root) | (g$E$nums[, 1] == root))
      MST <- g
      if (any(g$E$todel ))
        MST <- deleteEdges(g)
    } else {
      stop("vine not implemented")
    }
  } else {
    MST <- g

    # get edge list
    edgesList <- g$E$nums
    uid <- sort(unique(as.vector(g$E$nums)))
    luid <- length(uid)

    if (luid > 2) {
      # transform to adjacency list
      adjacencyList <- lapply(uid, function(u)
        c(edgesList[edgesList[,1] == u,2],
          edgesList[edgesList[,2] == u,1]))

      # find a tree by traversing the graph
      dfsorder <- dfs(adjacencyList, 1)
      newEdgesList <- t(apply(dfsorder$E, 1, sort))

      ## delete unused edges
      MST$E$todel <- !duplicated(rbind(newEdgesList,
                                       edgesList))[-seq(1,luid-1)]
    }

    if (any(MST$E$todel))
      MST <- deleteEdges(MST)

  }

  ## return result
  MST
}


## fit pair-copulas for the first vine tree
fit.FirstTreeCopulas <- function(MST,
                                 U,
                                 X,
                                 formula,
                                 familyset,
                                 selectioncrit,
                                 indeptest,
                                 trunclevel,
                                 control,
                                 light,
                                 cores = 1) {

  ## initialize estimation results with empty list
  d <- nrow(MST$E$nums)
  pc.data <- lapply(1:d, function(i) NULL)

  ## prepare for estimation and store names
  for (i in 1:d) {
    ## get edge and corresponding data
    a <- MST$E$nums[i, ]
    pc.data[[i]]$zr1 <- U[, a[1]]
    pc.data[[i]]$zr2 <- U[, a[2]]

    ## set names for this edge
    if (is.null(MST$V$names[a[1]])) {
      MST$E$Copula.CondName.1[i] <- a[1]
    } else {
      MST$E$Copula.CondName.1[i] <- MST$V$names[a[1]]
    }
    if (is.null(MST$V$names[a[2]])) {
      MST$E$Copula.CondName.2[i] <- a[2]
    } else {
      MST$E$Copula.CondName.2[i] <- MST$V$names[a[2]]
    }
    if (is.null(MST$V$names[a[1]]) || is.null(MST$V$names[a[2]])) {
      MST$E$Copula.Name[i] <- paste(a[1], a[2], sep = " , ")
    } else {
      MST$E$Copula.Name[i] <- paste(MST$V$names[a[1]],
                                    MST$V$names[a[2]],
                                    sep = " , ")
    }

  }

  if (!is.na(trunclevel) & trunclevel == 0) {
    indeptest <- -Inf
  }

  ## estimate parameters and select family
  pc.fits <- lapply(1:length(pc.data), function(k) fit.ACopula(u1 = pc.data[[k]]$zr1,
                                                               u2 = pc.data[[k]]$zr2,
                                                               X = X,
                                                               formula[k],
                                                               familyset,
                                                               selectioncrit,
                                                               indeptest,
                                                               control,
                                                               light,
                                                               cores))


  ## store estimated model and pseudo-observations for next tree
  for (i in 1:d) {

    MST$E$fits[[i]] <- pc.fits[[i]]
    MST$E$Copula.CondData.1[i] <- list(pc.fits[[i]]$CondOn.1)
    MST$E$Copula.CondData.2[i] <- list(pc.fits[[i]]$CondOn.2)


  }

  ## return results
  MST
}

## fit pair-copulas for vine trees > 1
fit.TreeCopulas <- function(X,
                            MST,
                            oldVineGraph,
                            formula,
                            familyset,
                            selectioncrit,
                            indeptest,
                            trunclevel,
                            control,
                            light,
                            cores = 1) {

  ## initialize estimation results with empty list
  d <- nrow(MST$E$nums)
  pc.data <- lapply(1:d, function(i) NULL)

  ## prepare for estimation
  for (i in 1:d) {
    ## get edge and corresponding data
    con <- MST$E$nums[i, ]
    temp <- oldVineGraph$E$nums[con, ]

    ## fetch corresponding data and names
    if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
      same <- temp[2, 1]
    } else {
      if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
        same <- temp[2, 2]
      }
    }

    if (temp[1, 1] == same) {
      zr1 <- oldVineGraph$E$Copula.CondData.2[con[1]]
      n1  <- oldVineGraph$E$Copula.CondName.2[con[1]]
    } else {
      zr1 <- oldVineGraph$E$Copula.CondData.1[con[1]]
      n1  <- oldVineGraph$E$Copula.CondName.1[con[1]]
    }
    if (temp[2, 1] == same) {
      zr2 <- oldVineGraph$E$Copula.CondData.2[con[2]]
      n2  <- oldVineGraph$E$Copula.CondName.2[con[2]]
    } else {
      zr2 <- oldVineGraph$E$Copula.CondData.1[con[2]]
      n2  <- oldVineGraph$E$Copula.CondName.1[con[2]]
    }

    if (is.list(zr1)) {
      zr1a <- as.vector(zr1[[1]])
      zr2a <- as.vector(zr2[[1]])
      n1a <- as.vector(n1[[1]])
      n2a <- as.vector(n2[[1]])
    } else {
      zr1a <- zr1
      zr2a <- zr2
      n1a <- n1
      n2a <- n2
    }

    pc.data[[i]]$zr1 <- zr1a
    pc.data[[i]]$zr2 <- zr2a

    MST$E$Copula.CondName.1[i] <- n1a
    MST$E$Copula.CondName.2[i] <- n2a
  }

  ## estimate parameters and select family
  pc.fits <- lapply(1:length(pc.data), function(k) fit.ACopula(u1 = pc.data[[k]]$zr1,
                                                               u2 = pc.data[[k]]$zr2,
                                                               X = X,
                                                               formula[k],
                                                               familyset,
                                                               selectioncrit,
                                                               indeptest,
                                                               control,
                                                               light,
                                                               cores))


  ## store estimated model and pseudo-observations for next tree
  for (i in 1:d) {

    MST$E$fits[[i]] <- pc.fits[[i]]
    MST$E$Copula.CondData.1[i] <- list(pc.fits[[i]]$CondOn.1)
    MST$E$Copula.CondData.2[i] <- list(pc.fits[[i]]$CondOn.2)

  }

  ## return results
  MST
}

## bivariate copula selection
fit.ACopula <- function(u1,
                        u2,
                        X,
                        formula,
                        familyset = NA,
                        selectioncrit = "aic",
                        indeptest = NA,
                        control = boost_control(),
                        light = FALSE,
                        cores = 1) {


  ## select family and estimate parameter(s) for the pair copula
  out <- suppressWarnings(boostBiCopSelect(formula = as.formula(formula),
                                           U = cbind(u1, u2),
                                           X = X,
                                           familyset = familyset,
                                           selectioncrit = selectioncrit,
                                           indeptest = indeptest,
                                           control = control,
                                           na.action = na.pass,
                                           light = light,
                                           cores = cores))

  ## change rotation if family is not symmetric wrt the main diagonal
  # parameter
  par1 <- unlist(predict(out, newdata = X, type = "parameter")$parameter)
  par2 <- out$par2
  # family
  family <- getFams(as.numeric(out$family))
  fam <- rep(family[1], length(par1))
  if (length(family) != 1) {
    fam[par1 < 0] <- family[2]
  }

  if (any(fam %in% c(23, 24, 26:30))) {
    fam[fam %in% c(23, 24, 26:30)] <- fam[fam %in% c(23, 24, 26:30)] + 10
  } else if (any(fam %in% c(33, 34, 36:40))) {
    fam[fam %in% c(33, 34, 36:40)] <- fam[fam %in% c(33, 34, 36:40)] - 10
  }


  ## store pseudo-observations for estimation in next tree
  out$CondOn.1 <- suppressWarnings(BiCopHfunc1(u1 = u2,
                                               u2 = u1,
                                               family = fam,
                                               par = par1,
                                               par2 = par2,
                                               check.pars = FALSE))
  out$CondOn.2 <- suppressWarnings(BiCopHfunc2(u1 = u2,
                                               u2 = u1,
                                               family = fam,
                                               par = par1,
                                               par2 = par2,
                                               check.pars = FALSE))
  ## return results
  return(out)
}

## build R-Vine matrix object based on nested set of trees
as.boostRVineCopula <- function(RVine, data, X, settings) {

  ## initialize objects
  n <- length(RVine$Tree) + 1
  nam <- RVine$Tree[[1]]$V$names
  nedSets <- list()
  crspfits <- list()

  ## get selected pairs, families and estimated parameters
  for (k in 1:(n - 2)) {

    nedSets[[k]]    <- RVine$Tree[[k]]$E$conditionedSet
    crspfits[[k]]   <- as.list(RVine$Tree[[k]]$E$fits)

  }
  crspfits[[n - 1]]   <- as.list(RVine$Tree[[n - 1]]$E$fits)
  if (is.list(RVine$Tree[[1]]$E$conditionedSet)) {
    nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet[[1]])
  } else {
    nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet)
  }

  ## initialize matrices for RVineMatrix object
  M <- matrix(NA, n, n)
  pair_copulas <- lapply(1:n, function(k) lapply(1:n, function(j) NULL))

  ## store structure, families and parameters in matrices
  for (k in 1:(n - 1)) {
    w <- nedSets[[n - k]][[1]][1]

    M[k, k] <- w
    M[(k + 1), k] <- nedSets[[n - k]][[1]][2]
    pair_copulas[[(k + 1)]][[k]]  <- crspfits[[n - k]][[1]]
    settings$boostRVM$family[k + 1, k] <- as.numeric(pair_copulas[[k + 1]][[k]]$family)

    if (k == (n - 1)) {
      M[(k + 1), (k + 1)] <- nedSets[[n - k]][[1]][2]
    } else {
      for (i in (k + 2):n) {
        for (j in 1:length(nedSets[[n - i + 1]])) {
          cs <- nedSets[[n - i + 1]][[j]]
          if (cs[1] == w) {
            M[i, k] <- cs[2]
            break
          } else if (cs[2] == w) {
            M[i, k] <- cs[1]
            break
          }
        }

        pair_copulas[[i]][[k]] <- crspfits[[n - i + 1]][[j]]
        settings$boostRVM$family[i, k] <- as.numeric(pair_copulas[[i]][[k]]$family)
        nedSets[[n - i + 1]][[j]]    <- NULL
        crspfits[[n - i + 1]][[j]] <- NULL

      }
    }
  }

  ## clean NAs
  M[is.na(M)] <- 0
  settings$boostRVM$Matrix <- M


  # calculate log-likelihood of model
  d <- ncol(settings$boostRVM$Matrix)
  l <- unlist(lapply(1:d, function(k) Filter(Negate(is.null), pair_copulas[[k]])), recursive = F)
  loglik <- sum(sapply(l, logLik))
  df <- sum(sapply(1:(d*(d-1)/2), function(k) {attr(logLik(l[[k]]), "df")}))
  aic <- sum(sapply(l, AIC))

  out <- list(model_frame = data.frame(data, X),
              pair_copulas = pair_copulas,
              boostRVM = settings$boostRVM,
              control = settings$control,
              stats = list(n = nrow(data),
                           df = df,
                           loglik = loglik,
                           aic = aic))

  class(out) <- "boostRVineCopula"

  return(out)
}

## functions for handling the tree structure
graphFromWeightMatrix <- function(W) {
  d <- ncol(W)
  # get variable names
  nms <- colnames(W)
  if (is.null(nms))
    nms <- paste0("V", 1:d)
  # construct edge set
  E <- cbind(do.call(c, sapply(1:(d-1), function(i) seq.int(i))),
             do.call(c, sapply(1:(d-1), function(i) rep(i+1, i))))
  # add edge names
  E.names <- apply(E, 1, function(x) paste(nms[x[1]],  nms[x[2]], sep = ","))
  # set weights
  w <- W[upper.tri(W)]

  ## return results
  list(V = list(names = nms,
                conditionedSet = NULL,
                conditioningSet = NULL),
       E = list(nums = E,
                names = E.names,
                weights = w,
                conditionedSet = lapply(1:nrow(E), function(i) E[i, ]),
                conditioningSet = NULL))
}
makeFullGraph <- function(d) {
  ## create matrix of all combinations
  E <- cbind(do.call(c, lapply(1:(d-1), function(i) rep(i, d-i))),
             do.call(c, lapply(1:(d-1), function(i) (i+1):d)))
  E <- matrix(E, ncol = 2)

  ## output dummy list with edges set
  list(V = list(names = NULL,
                conditionedSet = NULL,
                conditioningSet = NULL),
       E = list(nums = E,
                names = NULL,
                weights = NULL,
                conditionedSet = E,
                conditioningSet = NULL))
}
adjacencyMatrix <- function(g) {
  ## create matrix of all combinations
  d <- length(g$V$names)
  v.all <- cbind(do.call(c, lapply(1:(d-1), function(i) seq.int(i))),
                 do.call(c, lapply(1:(d-1), function(i) rep(i+1, i))))

  ## find weight
  vals <- apply(v.all, 1, set_weight, E = g$E)

  ## create symmetric matrix of weights
  M <- matrix(0, d, d)
  M[upper.tri(M)] <- vals
  M <- M + t(M)
  diag(M) <- Inf

  ## return final matrix
  M
}
set_weight <- function(x, E) {
  ## convert weights so that minimum spanning tree algorithm can be applied
  is.edge <- (x[1] == E$nums[, 1]) & (x[2] == E$nums[, 2])
  if (!any(is.edge)) Inf else -E$weights[which(is.edge)]
}
deleteEdges <- function(g) {
  ## reduce edge list
  keep <- which(!g$E$todel)
  E <- list(nums            = matrix(g$E$nums[keep, ], ncol = 2),
            names           = g$E$names[keep],
            weights         = g$E$weights[keep],
            conditionedSet  = g$E$conditionedSet[keep],
            conditioningSet = g$E$conditioningSet[keep])

  ## return reduced graph
  list(V = g$V, E = E)
}
## initialize graph for next vine tree (possible edges)
buildNextGraph <- function(oldVineGraph, treecrit, truncated = FALSE) {

  d <- nrow(oldVineGraph$E$nums)

  ## initialize with full graph
  g <- makeFullGraph(d)
  g$V$names <- oldVineGraph$E$names
  g$V$conditionedSet <- oldVineGraph$E$conditionedSet
  g$V$conditioningSet <- oldVineGraph$E$conditioningSet

  ## get info for all edges
  out <- lapply(seq_len(nrow(g$E$nums)),
                getEdgeInfo,
                g = g,
                oldVineGraph = oldVineGraph,
                treecrit = treecrit,
                truncated = truncated)

  ## annotate graph (same order as in old version of this function)
  g$E$weights         <- sapply(out, function(x) x$w)
  g$E$names           <- sapply(out, function(x) x$name)
  g$E$conditionedSet  <- lapply(out, function(x) x$nedSet)
  g$E$conditioningSet <- lapply(out, function(x) x$ningSet)
  g$E$todel           <- sapply(out, function(x) x$todel)

  ## delete edges that are prohibited by the proximity condition
  deleteEdges(g)
}
## function for obtaining edge information
getEdgeInfo <- function(i, g, oldVineGraph, treecrit, truncated = FALSE) {

  ## get edge
  con <- g$E$nums[i, ]
  temp <- oldVineGraph$E$nums[con, ]


  ## check for proximity condition
  ok <- FALSE
  if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
    ok <- TRUE
    same <- temp[2, 1]
  } else {
    if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
      ok <- TRUE
      same <- temp[2, 2]
    }
  }

  ## dummy output
  w <- nedSet <- ningSet <- name <- NA
  todel <- TRUE

  # info if proximity condition is fulfilled ...
  if (ok) {
    ## infer conditioned set and conditioning set
    l1 <- c(g$V$conditionedSet[[con[1]]],
            g$V$conditioningSet[[con[1]]])
    l2 <- c(g$V$conditionedSet[[con[2]]],
            g$V$conditioningSet[[con[2]]])
    nedSet <- c(setdiff(l1, l2), setdiff(l2, l1))
    ningSet <- intersect(l1, l2)

    ## mark as ok
    todel <- FALSE

    if (truncated == FALSE) {
      ## get data
      if (temp[1, 1] == same) {
        zr1 <- oldVineGraph$E$Copula.CondData.2[con[1]]
      } else {
        zr1 <- oldVineGraph$E$Copula.CondData.1[con[1]]
      }
      if (temp[2, 1] == same) {
        zr2 <- oldVineGraph$E$Copula.CondData.2[con[2]]
      } else {
        zr2 <- oldVineGraph$E$Copula.CondData.1[con[2]]
      }
      if (is.list(zr1)) {
        zr1a <- as.vector(zr1[[1]])
        zr2a <- as.vector(zr2[[1]])
      } else {
        zr1a <- zr1
        zr2a <- zr2
      }

      ## calculate Kendall's tau
      keine_nas <- !(is.na(zr1a) | is.na(zr2a))
      w <- treecrit(zr1a[keine_nas], zr2a[keine_nas])

      ## get names
      name.node1 <- strsplit(g$V$names[con[1]], split = " *[,;] *")[[1]]
      name.node2 <- strsplit(g$V$names[con[2]], split = " *[,;] *")[[1]]

      ## set edge name
      nmdiff <- c(setdiff(name.node1, name.node2),
                  setdiff(name.node2, name.node1))
      nmsect <- intersect(name.node1, name.node2)
      name <- paste(paste(nmdiff, collapse = ","),
                    paste(nmsect, collapse = ","),
                    sep = " ; ")
    } else {
      w <- 1
    }
  }

  ## return edge information
  list(w = w,
       nedSet = nedSet,
       ningSet = ningSet,
       name = name,
       todel = todel)
}

# depth first search to build a tree without finding the MST
dfs <- function(adjacencyList, v, e = NULL, dfsorder = list()) {
  dfsorder$V <- c(dfsorder$V, v)
  dfsorder$E <- rbind(dfsorder$E, e)
  for (u in adjacencyList[[v]]) {
    if (!is.element(u, dfsorder$V)) {
      dfsorder <- dfs(adjacencyList, u, c(u,v), dfsorder)
    }
  }
  return(dfsorder)
}


