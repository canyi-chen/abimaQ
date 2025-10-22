#' @title PDF/PMF for General GLM Families
#' @description
#' Computes \eqn{f_Y(y \mid x)} for standard GLM families (continuous densities
#' or discrete mass functions). Uses native formulas or SPA for continuous EDFs.
#'
#' @inheritParams cdf_glm_general
#'
#' @return Numeric vector of PDF or PMF values.
#' @examples
#' x <- seq(0, 10, by = 0.5)
#' pdf_glm_general(x, eta = 2, family = Gamma(link = "log"), dispersion = 1)
#' @export
pdf_glm_general <- function(
    y, eta, family,
    dispersion = 1,
    method = c("auto","native","saddlepoint"),
    nb_size = NULL
) {
  method <- match.arg(method)
  y   <- as.numeric(y)
  eta <- as.numeric(eta)
  if (length(eta) == 1L && length(y) > 1L) eta <- rep(eta, length(y))
  if (length(y) != length(eta)) stop("y and eta length mismatch.")

  fam  <- family$family
  link <- family$link

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

  pdf_native <- function(y, mu) {
    switch(fam,
           "gaussian" = stats::dnorm(y, mean = mu, sd = sqrt(dispersion)),
           "binomial" = stats::dbinom(y, size = 1L, prob = pmin(pmax(mu,1e-12),1-1e-12)),
           "poisson"  = stats::dpois(floor(y), lambda = mu),
           "Gamma"    = { shape <- 1/dispersion; scale <- dispersion * mu
           stats::dgamma(y, shape = shape, scale = scale) },
           "inverse.gaussian" = {
             lambda <- 1/dispersion
             dens <- rep(0, length(y)); ok <- (y > 0)
             dens[ok] <- sqrt(lambda/(2*pi*y[ok]^3)) *
               exp(-lambda*(y[ok]-mu[ok])^2 / (2*mu[ok]^2*y[ok]))
             dens
           },
           {
             size <- nb_size; if (is.null(size)) stop("Provide nb_size for NB density.")
             stats::dnbinom(floor(y), size = size, mu = mu)
           }
    )
  }

  if (method %in% c("auto","native") && native_ok()) return(pdf_native(y, mu))

  # SPA density for continuous EDFs
  fam_specs <- switch(fam,
                      "gaussian" = list(
                        theta = function(mu) mu, b = function(th) 0.5*th^2,
                        b1 = function(th) th,   b2 = function(th) 1
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
  if (is.null(fam_specs))
    stop("Saddlepoint density not implemented for this family; use method='native'.")

  aphi   <- dispersion
  theta0 <- fam_specs$theta(mu)
  K   <- function(t, i) (fam_specs$b(theta0[i]+t) - fam_specs$b(theta0[i]))/aphi
  K1  <- function(t, i) (fam_specs$b1(theta0[i]+t))/aphi
  K2  <- function(t, i) (fam_specs$b2(theta0[i]+t))/aphi

  t_hat <- numeric(length(y))
  for (i in seq_along(y)) {
    target <- y[i]; f <- function(t) K1(t,i) - target
    lo <- -5; hi <- 5
    for (k in 1:6) { if (f(lo)*f(hi) < 0) break; lo <- lo*2; hi <- hi*2 }
    t_hat[i] <- tryCatch(uniroot(f, lo, hi, tol = 1e-10)$root, error = function(e) NA_real_)
  }

  dens <- rep(NA_real_, length(y))
  for (i in seq_along(y)) {
    th <- t_hat[i]; if (!is.finite(th)) next
    K2s <- K2(th,i); if (!is.finite(K2s) || K2s <= 0) next
    dens[i] <- exp(K(th,i) - th*y[i]) / sqrt(2*pi*K2s)
  }
  miss <- !is.finite(dens)
  if (any(miss)) dens[miss] <- pdf_native(y[miss], mu[miss])
  dens
}
