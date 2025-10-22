#' CDF for General GLM Families
#'
#' Computes \eqn{F_Y(y \mid x)} for common GLM families (Gaussian, Binomial,
#' Poisson, Gamma, Inverse-Gaussian, Negative-Binomial) using native formulas
#' or a saddlepoint approximation (SPA), with a small MC fallback.
#'
#' @param y Numeric vector of evaluation points.
#' @param eta Numeric vector of linear predictors \eqn{x^\top \beta}.
#' @param family A \code{\link[stats]{family}} object.
#' @param dispersion Numeric dispersion \eqn{\phi}.
#' @param method One of \code{"auto"}, \code{"native"}, \code{"saddlepoint"}.
#' @param discrete_pit For discrete families, one of \code{"midp"}, \code{"jitter"}, \code{"closed"}.
#' @param nb_size Optional Negative-Binomial size (theta).
#'
#' @return Numeric vector of CDF values.
#' @examples
#' x <- seq(0, 10, by = 0.5)
#' cdf_glm_general(x, eta = 2, family = Gamma(link = "log"), dispersion = 1)
#' @export
cdf_glm_general <- function(
    y, eta, family,
    dispersion = 1,
    method = c("auto","native","saddlepoint"),
    discrete_pit = c("midp","jitter","closed"),
    nb_size = NULL
) {
  method <- match.arg(method)
  discrete_pit <- match.arg(discrete_pit)
  y   <- as.numeric(y)
  eta <- as.numeric(eta)
  if (length(eta) == 1L && length(y) > 1L) eta <- rep(eta, length(y))
  if (length(y) != length(eta)) stop("y and eta must match in length (or eta length 1).")

  fam  <- family$family
  link <- family$link
  eps  <- 1e-12
  clip01 <- function(p) pmin(pmax(p, eps), 1 - eps)

  invlink <- switch(link,
                    "identity" = function(e) e,
                    "log"      = function(e) exp(e),
                    "inverse"  = function(e) 1/e,
                    "logit"    = function(e) 1/(1+exp(-e)),
                    "probit"   = function(e) stats::pnorm(e),
                    "cloglog"  = function(e) 1 - exp(-exp(e)),
                    function(e) exp(e)
  )
  mu <- invlink(eta)

  native_ok <- function() {
    switch(fam,
           "gaussian" = TRUE,
           "binomial" = TRUE,
           "poisson"  = TRUE,
           "Gamma"    = TRUE,
           "inverse.gaussian" = TRUE,
           { grepl("Negative Binomial", fam, ignore.case = TRUE) }
    )
  }

  # Inverse-Gaussian CDF
  pinvgauss_cf <- function(yy, mu, lambda) {
    z1 <- sqrt(lambda/yy) * (yy/mu - 1)
    z2 <- -sqrt(lambda/yy) * (yy/mu + 1)
    stats::pnorm(z1) + exp(2*lambda/mu) * stats::pnorm(z2)
  }

  cdf_native <- function(y, mu) {
    switch(fam,
           "gaussian" = stats::pnorm(y, mean = mu, sd = sqrt(dispersion)),
           "binomial" = {
             p1 <- clip01(mu)
             Fy <- ifelse(y < 1, 1 - p1, 1)
             if (discrete_pit == "closed") return(Fy)
             Fy_minus <- ifelse(y < 1, 0, 1 - p1)
             if (discrete_pit == "midp") return(clip01(0.5 * (Fy_minus + Fy)))
             Fy_minus + (Fy - Fy_minus) * stats::runif(length(y))
           },
           "poisson" = {
             Fy <- stats::ppois(floor(y), lambda = mu)
             if (discrete_pit == "closed") return(Fy)
             Fy_minus <- stats::ppois(pmax(floor(y) - 1, -1), lambda = mu)
             if (discrete_pit == "midp") return(clip01(0.5 * (Fy_minus + Fy)))
             Fy_minus + (Fy - Fy_minus) * stats::runif(length(y))
           },
           "Gamma" = {
             shape <- 1/dispersion; scale <- dispersion * mu
             stats::pgamma(y, shape = shape, scale = scale)
           },
           "inverse.gaussian" = {
             lambda <- 1/dispersion
             pinvgauss_cf(y, mu = mu, lambda = lambda)
           },
           {
             size <- nb_size
             if (is.null(size)) stop("Provide nb_size for Negative Binomial CDF.")
             stats::pnbinom(floor(y), size = size, mu = mu)
           }
    )
  }

  if (method %in% c("auto","native") && native_ok()) return(cdf_native(y, mu))

  # Lugannani–Rice SPA (continuous EDFs)
  fam_specs <- switch(fam,
                      "gaussian" = list(
                        theta = function(mu) mu, b = function(th) 0.5*th^2,
                        b1 = function(th) th, b2 = function(th) 1
                      ),
                      "poisson" = list(
                        theta = function(mu) log(mu), b = function(th) exp(th),
                        b1 = function(th) exp(th),    b2 = function(th) exp(th)
                      ),
                      "binomial" = list(
                        theta = function(mu) log(mu/(1-mu)),
                        b  = function(th) log1p(exp(th)),
                        b1 = function(th) 1/(1+exp(-th)),
                        b2 = function(th) exp(th)/(1+exp(th))^2
                      ),
                      "Gamma" = list(
                        theta = function(mu) -1/mu, b = function(th) -log(-th),
                        b1 = function(th) -1/th,    b2 = function(th) 1/th^2
                      ),
                      "inverse.gaussian" = list(
                        theta = function(mu) -1/(2*mu^2), b = function(th) -sqrt(-2*th),
                        b1 = function(th) -1/sqrt(-2*th),
                        b2 = function(th) -1/(2*(-2*th)^(3/2)) * (-2)
                      ),
                      NULL
  )
  if (is.null(fam_specs)) return(cdf_mc(y, mu, fam))  # last resort

  aphi   <- dispersion
  theta0 <- fam_specs$theta(mu)
  K   <- function(t, i) (fam_specs$b(theta0[i] + t) - fam_specs$b(theta0[i])) / aphi
  K1  <- function(t, i) (fam_specs$b1(theta0[i] + t)) / aphi
  K2  <- function(t, i) (fam_specs$b2(theta0[i] + t)) / aphi

  t_hat <- numeric(length(y))
  for (i in seq_along(y)) {
    target <- y[i]
    f  <- function(t) K1(t, i) - target
    lo <- -5; hi <- 5
    for (k in 1:6) { if (f(lo)*f(hi) < 0) break; lo <- lo*2; hi <- hi*2 }
    t_hat[i] <- tryCatch(uniroot(f, lo, hi, tol = 1e-10)$root, error = function(e) NA_real_)
  }

  w <- u <- rep(NA_real_, length(y))
  for (i in seq_along(y)) {
    th <- t_hat[i]
    Ks  <- K(th, i)
    K2s <- K2(th, i)
    if (!is.finite(th) || !is.finite(Ks) || !is.finite(K2s) || K2s <= 0) next
    w[i] <- sign(th) * sqrt( max(2*(th*y[i] - Ks), 0) )
    u[i] <- th * sqrt(K2s)
  }

  phiN <- function(z) exp(-0.5*z^2)/sqrt(2*pi)
  Fspa <- stats::pnorm(w) + phiN(w) * (1/w - 1/u)
  Fspa[!is.finite(Fspa)] <- NA_real_
  miss <- is.na(Fspa)
  if (any(miss)) Fspa[miss] <- cdf_mc(y[miss], mu[miss], fam)
  pmin(pmax(Fspa, eps), 1 - eps)
}
