#' Plot Coefficient Paths of boostBiCop object.
#'
#' @param x an x of class \code{boostBiCop}.
#' @param col character; colour(s) for the coefficient paths.
#' @param loglik logical; either \code{FALSE} (default), then the coefficient paths are plotted.
#' Otherwise, if \code{TRUE}, then the log-likelihood contribution paths are plotted.
#' @param ... unused.
#'
#' @importFrom graphics axis par strwidth grid
#' @importFrom stats ts
#' @importFrom colorspace qualitative_hcl
#'
#' @examples
#' # see boostBiCopEst
#'
#' @export
plot.boostBiCop <- function(x, col, loglik = FALSE, ...) {

  if (length(x$paths$cov_path) == 1) {
    stop("Number of boosting iterations is to small!")
  }

  if (!loglik) {

    path <- x$paths$coef_path
    x_names <- colnames(path)
    x_sel <- sort(unique(x$paths$cov_path[-1]))
    path <- as.matrix(path[, x_names[x_sel]])
    colnames(path) <- x_names[x_sel]

    # plot settings
    main <- "Coefficient path"
    xlab <- "Boosting iterations"
    ylab <- "Coefficient"

  } else {

    x_names <- names(x$coefficients)
    loglik_diff <- diff(x$paths$loglik_path)
    cov_path <- x$paths$cov_path[-1]
    covs <- sort(unique(cov_path))
    path <- matrix(0, nrow = x$iterations+1, ncol = length(x$coefficients), dimnames = list(c(), names(x$coefficients)))
    for (k in 2:(x$iterations+1)) {
      path[k, cov_path[k-1]] <- path[k-1, cov_path[k-1]] + loglik_diff[k-1]
      path[k, -cov_path[k-1]] <- path[k-1, -cov_path[k-1]]
    }
    if (any(is.na(cov_path))) {
      index_NA <- which(is.na(cov_path))
      path[index_NA+1, ] <- matrix(rep(path[index_NA[1], ], length(index_NA)), ncol = length(x$coefficients), byrow = T)
    }
    path <- as.matrix(path[, covs])
    colnames(path) <- x_names[covs]

    # plot settings
    main <- "Coefficient path"
    xlab <- "Boosting iterations"
    ylab <- "log-likelihood contribution"

  }

  if (missing(col)) {
    # to 'achieve' a 'better' color diversity in the visualization (especially in high-dim. settings)
    set.seed(123)
    col <- sample(qualitative_hcl(ncol(path), palette = "Set 2"), size = ncol(path), replace = F)
  }

  # adjust figure margins
  par(mai = c(1, 1, 1, max(strwidth(colnames(path), units = "in")) + 0.5))

  # plot coefficient path
  plot(ts(path, start = 0),
       plot.type = "single",
       main = main,
       xlab = xlab,
       ylab = ylab,
       col = col,
       type = "s")
  axis(4,
       at = as.numeric(path[x$iterations+1, ]),
       labels = colnames(path),
       las = 1,
       col.axis = 1,
       col.ticks = 1)
  grid()


}



