% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boostBiCopSim.R
\name{boostBiCopSim}
\alias{boostBiCopSim}
\title{Simulation from a Conditional Bivariate Copula}
\usage{
boostBiCopSim(object, N, X, ...)
}
\arguments{
\item{object}{an \code{boostBiCop} object.}

\item{N}{integer, specifying the number of simulation samples.}

\item{X}{matrix containing the covariates.}

\item{...}{unused.}
}
\value{
A matrix with two columns containing the samples of the \code{boostBiCop} object, evaluated at \code{X}.
}
\description{
Simulation from a \code{boostBiCop} object.
}
\examples{
# load simulated data
data(data_bicop)

# fit object
object <- boostBiCopEst(formula = ~.,
                        U = data_bicop[, 1:2],
                        X = data_bicop[, -c(1:2)],
                        family = 301,
                        control = boost_control(deselection = "attributable"))
boostBiCopSim(object,
              N = 10,
              X = data_bicop[, -c(1:2)])

}
