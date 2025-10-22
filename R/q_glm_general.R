#' @title Quantile Function for GLM Families
#' @description
#' Computes \eqn{F_Y^{-1}(p \mid x)}. Uses native quantiles when available,
#' otherwise numerical inversion of the CDF.
#'
#' @param p Numeric vector of probabilities in \eqn{(0,1)}.
#' @inheritParams cdf_glm_general
#' @param lower,upper Bounds for root finding.
#'
#' @return Numeric vector of quantiles.
#' @examples
#' q_glm_general(0.5, eta = 2, family = Gamma(link = "log"), dispersion = 1)
#' @export
q_glm_general <- function(
    p, eta, family, dispersion = 1, nb_size = NULL,
    lower = .Machine$double.eps, upper = 1e6
) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  fam <- family$family
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

  if (fam == "gaussian")  return(stats::qnorm(p, mean = mu, sd = sqrt(dispersion)))
  if (fam == "poisson")   return(stats::qpois(p, lambda = mu))
  if (fam == "binomial")  return(stats::qbinom(p, size = 1L, prob = pmin(pmax(mu,1e-12),1-1e-12)))
  if (fam == "Gamma")     { shape <- 1/dispersion; scale <- dispersion * mu; return(stats::qgamma(p, shape=shape, scale=scale)) }
  if (grepl("Negative Binomial", fam, ignore.case = TRUE)) {
    size <- nb_size; if (is.null(size)) stop("Provide nb_size for Negative Binomial quantile.")
    return(stats::qnbinom(p, size = size, mu = mu))
  }
  f <- function(q) cdf_glm_general(q, eta, family, dispersion, method="auto", nb_size=nb_size) - p
  uniroot(f, lower = lower, upper = upper, tol = 1e-10)$root
}
