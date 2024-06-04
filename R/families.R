## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib boostCopula, .registration = TRUE
## usethis namespace: end
NULL

#' Copula families and corresponding helper function
#'
#' @importFrom VineCopula BiCopPDF BiCopDeriv BiCopEst
#' @noRd

##################
# helper functions
##################
tau2par <- function(family, tau) {

  if (family == 0) {
    par <- rep(0, length(tau))
  }
  if (family == 1) {
    # correct to avoid problems with abs(tau) = 1
    tau <- pmin(pmax(tau, -1 + 1e-3), 1 - 1e-3)
    par <- sin(pi/2*tau)
  }
  if (family %in% c(301:304)) {
    # correct to avoid problems with abs(tau) = 1
    tau <- pmin(pmax(tau, -1 + 1e-7), 1 - 1e-7)
    par <- 2*tau/(1 - abs(tau))
    # parameter needs to be greater 0
    par[par == 0] <- 1e-7
  }
  if (family %in% c(401:404)) {
    # correct to avoid problems with abs(tau) = 1
    tau <- pmin(pmax(tau, -1 + 1e-7), 1 - 1e-7)
    par <- sign(tau)/((1 - abs(tau)))
    # correct if sign(tau) = 0
    par[par == 0] <- 1
  }
  return(par)
}
getFams <- function(family) {
  if (family == 0) {
    fam <- 0
  } else if (family == 1) {
    fam <- 1
  } else if (family %in% c(301:304)) {
    fam <- as.numeric(rev(expand.grid(c(23, 33), c(3, 13)))[family - 300, ])
  } else if (family %in% c(401:404)) {
    fam <- as.numeric(rev(expand.grid(c(24, 34), c(4, 14)))[family - 400, ])
  }
  return(fam)
}

#################
# Copula families
#################
BiCop_1_tau <- function() {
  list(
    loss = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-3), 1 - 1e-3)
      par <- sin(pi/2 * f)

      out <- - 1*log(BiCopPDF(u1 = u1,
                              u2 = u2,
                              family = 1,
                              par = par,
                              par2 = par2,
                              check.pars = FALSE))
      out[is.nan(out)] <- Inf
      return(out)

    },
    ngradient = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-3), 1 - 1e-3)
      par <- sin(pi/2 * f)

      out <- BiCopDeriv(u1 = u1,
                        u2 = u2,
                        family = 1,
                        par = par,
                        par2 = par2,
                        deriv = "par",
                        log = TRUE,
                        check.pars = FALSE) * 1/(2*pi)*sqrt(1 - par^2) * (pi^2 - 4*asin(par)^2)

      return(out)

    },
    family = "1",
    familyname = "Gaussian",
    par2 = 0
  )

}
BiCop_301_tau <- function(){

  list(
    loss = function(u, f, par2){

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- 2*f/(1 - abs(f))
      par[par == 0] <- 1e-7
      fam <- rep(3, length(par))
      fam[par < 0] <- 23

      out <- -1*log(BiCopPDF(u1 = u1,
                             u2 = u2,
                             family = fam,
                             par = par,
                             par2 = par2,
                             check.pars = FALSE))

      out[is.nan(out)] <- Inf
      return(out)

    },
    ngradient = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- 2*f/(1 - abs(f))
      par[par == 0] <- 1e-7
      fam <- rep(3, length(par))
      fam[par < 0] <- 23

      out <- BiCopDeriv(u1 = u1,
                        u2 = u2,
                        family = fam,
                        par = par,
                        par2 = par2,
                        deriv = "par",
                        log = TRUE,
                        check.pars = FALSE) * (2*sign(par)*par + 2)

      return(out)


    },
    family = "301",
    familyname = "double Clayton I",
    par2 = 0
  )

}
BiCop_302_tau <- function(){

  list(
    loss = function(u, f, par2){

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- 2*f/(1 - abs(f))
      par[par == 0] <- 1e-7
      fam <- rep(3, length(par))
      fam[par < 0] <- 33

      out <- -1*log(BiCopPDF(u1 = u1,
                             u2 = u2,
                             family = fam,
                             par = par,
                             par2 = par2,
                             check.pars = FALSE))
      out[is.nan(out)] <- Inf
      return(out)

    },
    ngradient = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- 2*f/(1 - abs(f))
      par[par == 0] <- 1e-7
      fam <- rep(3, length(par))
      fam[par < 0] <- 33

      out <- BiCopDeriv(u1 = u1,
                        u2 = u2,
                        family = fam,
                        par = par,
                        par2 = par2,
                        deriv = "par",
                        log = TRUE,
                        check.pars = FALSE) * (2*sign(par)*par + 2)

      return(out)


    },
    family = "302",
    familyname = "double Clayton II",
    par2 = 0
  )

}
BiCop_303_tau <- function(){

  list(
    loss = function(u, f, par2){

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- 2*f/(1 - abs(f))
      par[par == 0] <- 1e-7
      fam <- rep(13, length(par))
      fam[par < 0] <- 23

      out <- -1*log(BiCopPDF(u1 = u1,
                             u2 = u2,
                             family = fam,
                             par = par,
                             par2 = par2,
                             check.pars = FALSE))
      out[is.nan(out)] <- Inf
      return(out)

    },
    ngradient = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- 2*f/(1 - abs(f))
      par[par == 0] <- 1e-7
      fam <- rep(13, length(par))
      fam[par < 0] <- 23

      out <- BiCopDeriv(u1 = u1,
                        u2 = u2,
                        family = fam,
                        par = par,
                        par2 = par2,
                        deriv = "par",
                        log = TRUE,
                        check.pars = FALSE) * (2*sign(par)*par + 2)


      return(out)


    },
    family = "303",
    familyname = "double Clayton III",
    par2 = 0
  )

}
BiCop_304_tau <- function(){

  list(
    loss = function(u, f, par2){

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- 2*f/(1 - abs(f))
      par[par == 0] <- 1e-7
      fam <- rep(13, length(par))
      fam[par < 0] <- 33

      out <- -1*log(BiCopPDF(u1 = u1,
                             u2 = u2,
                             family = fam,
                             par = par,
                             par2 = par2,
                             check.pars = FALSE))
      out[is.nan(out)] <- Inf
      return(out)

    },
    ngradient = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- 2*f/(1 - abs(f))
      par[par == 0] <- 1e-7
      fam <- rep(13, length(par))
      fam[par < 0] <- 33

      out <- BiCopDeriv(u1 = u1,
                        u2 = u2,
                        family = fam,
                        par = par,
                        par2 = par2,
                        deriv = "par",
                        log = TRUE,
                        check.pars = FALSE) * (2*sign(par)*par + 2)

      return(out)


    },
    family = "304",
    familyname = "double Clayton IV",
    par2 = 0
  )

}
BiCop_401_tau <- function(){

  list(
    loss = function(u, f, par2){

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- sign(f)/((1 - abs(f)))
      # correct if sign(f) = 0
      par[par == 0] <- 1
      fam <- rep(4, length(par))
      fam[par < 0] <- 24

      out <- -1*log(BiCopPDF(u1 = u1,
                             u2 = u2,
                             family = fam,
                             par = par,
                             par2 = par2,
                             check.pars = FALSE))
      out[is.nan(out)] <- Inf
      return(out)

    },
    ngradient = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- sign(f)/((1 - abs(f)))
      # correct if sign(f) = 0
      par[par == 0] <- 1
      fam <- rep(4, length(par))
      fam[par < 0] <- 24

      out <- BiCopDeriv(u1 = u1,
                        u2 = u2,
                        family = fam,
                        par = par,
                        par2 = par2,
                        deriv = "par",
                        log = TRUE,
                        check.pars = FALSE) * (2*sign(par)*par - 1)

      return(out)


    },
    family = "401",
    familyname = "double Gumbel I",
    par2 = 0
  )

}
BiCop_402_tau <- function(){

  list(
    loss = function(u, f, par2){

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- sign(f)/((1 - abs(f)))
      # correct if sign(f) = 0
      par[par == 0] <- 1
      fam <- rep(4, length(par))
      fam[par < 0] <- 34

      out <- -1*log(BiCopPDF(u1 = u1,
                             u2 = u2,
                             family = fam,
                             par = par,
                             par2 = par2,
                             check.pars = FALSE))
      out[is.nan(out)] <- Inf
      return(out)

    },
    ngradient = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- sign(f)/((1 - abs(f)))
      # correct if sign(f) = 0
      par[par == 0] <- 1
      fam <- rep(4, length(par))
      fam[par < 0] <- 34

      out <- BiCopDeriv(u1 = u1,
                        u2 = u2,
                        family = fam,
                        par = par,
                        par2 = par2,
                        deriv = "par",
                        log = TRUE,
                        check.pars = FALSE) * (2*sign(par)*par - 1)

      return(out)


    },
    family = "402",
    familyname = "double Gumbel II",
    par2 = 0
  )

}
BiCop_403_tau <- function(){

  list(
    loss = function(u, f, par2){

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- sign(f)/((1 - abs(f)))
      # correct if sign(f) = 0
      par[par == 0] <- 1
      fam <- rep(14, length(par))
      fam[par < 0] <- 24

      out <- -1*log(BiCopPDF(u1 = u1,
                             u2 = u2,
                             family = fam,
                             par = par,
                             par2 = par2,
                             check.pars = FALSE))
      out[is.nan(out)] <- Inf
      return(out)

    },
    ngradient = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- sign(f)/((1 - abs(f)))
      # correct if sign(f) = 0
      par[par == 0] <- 1
      fam <- rep(14, length(par))
      fam[par < 0] <- 24

      out <- BiCopDeriv(u1 = u1,
                        u2 = u2,
                        family = fam,
                        par = par,
                        par2 = par2,
                        deriv = "par",
                        log = TRUE,
                        check.pars = FALSE) * (2*sign(par)*par - 1)

      return(out)


    },
    family = "403",
    familyname = "double Gumbel III",
    par2 = 0
  )

}
BiCop_404_tau <- function(){

  list(
    loss = function(u, f, par2){

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- sign(f)/((1 - abs(f)))
      # correct if sign(f) = 0
      par[par == 0] <- 1
      fam <- rep(14, length(par))
      fam[par < 0] <- 34

      out <- -1*log(BiCopPDF(u1 = u1,
                             u2 = u2,
                             family = fam,
                             par = par,
                             par2 = par2,
                             check.pars = FALSE))
      out[is.nan(out)] <- Inf
      return(out)

    },
    ngradient = function(u, f, par2) {

      u1 <- u[, 1]
      u2 <- u[, 2]

      f <- tanh(f)
      f <- pmin(pmax(f, -1 + 1e-7), 1 - 1e-7)
      par <- sign(f)/((1 - abs(f)))
      # correct if sign(f) = 0
      par[par == 0] <- 1
      fam <- rep(14, length(par))
      fam[par < 0] <- 34

      out <- BiCopDeriv(u1 = u1,
                        u2 = u2,
                        family = fam,
                        par = par,
                        par2 = par2,
                        deriv = "par",
                        log = TRUE,
                        check.pars = FALSE) * (2*sign(par)*par - 1)

      return(out)


    },
    family = "404",
    familyname = "double Gumbel IV",
    par2 = 0
  )

}
