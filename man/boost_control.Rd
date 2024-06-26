% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boost_control.R
\name{boost_control}
\alias{boost_control}
\title{Control Parameters for the Gradient-Boosting Algorithm}
\usage{
boost_control(
  maxit = 500,
  nu = 0.1,
  mstop = "aic",
  nfolds = 10,
  foldid = NULL,
  deselection = "none",
  gamma = 0.01,
  center = TRUE,
  cores.control = 1
)
}
\arguments{
\item{maxit}{integer; giving the number of initial boosting iterations. Default: \code{maxit = 500}.}

\item{nu}{double; defining the step size or shrinkage parameter. Default: \code{nu = 0.1}.}

\item{mstop}{character; stopping criterion "\code{max}", "\code{aic}" (default), "\code{cv}".}

\item{nfolds}{integer; number of folds in cross validation \code{mstop = "cv"}.}

\item{foldid}{vector; an optional vector of values between 1 and \code{nfolds} identifying the fold for each observation.
Default: \code{foldid = NULL}.}

\item{deselection}{character; "\code{attributable}", "\code{cumulative}" or "\code{none}" (default).}

\item{gamma}{double; defining the threshold for \code{deselection}. Default: \code{gamma = 0.01}.}

\item{center}{logical; indicating if the covariates should be mean centered \code{TRUE} (default) before fitting or not \code{FALSE}.}

\item{cores.control}{integer; the number of cores used for computations. Default: \code{cores = 1}.}
}
\value{
List with specified parameters.
}
\description{
Control Parameters for the Gradient-Boosting Algorithm
}
