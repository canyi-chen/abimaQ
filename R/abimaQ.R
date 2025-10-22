# =============================================================================
# quantile_mediation_glm.R  (package-ready)
# =============================================================================

# -------------------------------
# Variance & dispersion utilities
# -------------------------------

#' Variance Function V(mu) for Common GLM Families
#'
#' Returns the mean–variance relationship \eqn{V(\mu)} used by 1-parameter
#' exponential-family GLMs. Implemented families:
#' \itemize{
#'   \item \code{"gaussian"}: \eqn{V(\mu)=1}
#'   \item \code{"poisson"}: \eqn{V(\mu)=\mu}
#'   \item \code{"binomial"} (Bernoulli): \eqn{V(\mu)=\mu(1-\mu)}
#'   \item \code{"Gamma"}: \eqn{V(\mu)=\mu^2}
#'   \item \code{"inverse.gaussian"}: \eqn{V(\mu)=\mu^3}
#' }
#'
#' @param mu Numeric vector of means.
#' @param family_name Character scalar; one of the supported family names
#'   (typically \code{family$family} from a GLM fit).
#'
#' @return Numeric vector \code{V(mu)} of the same length as \code{mu}.
#'
#' @examples
#' .V_mu(mu = c(1,2,3), family_name = "poisson")
#'
#' @keywords internal
#' @noRd
.V_mu <- function(mu, family_name) {
  switch(family_name,
         "gaussian"          = rep(1, length(mu)),
         "poisson"           = mu,
         "binomial"          = mu * (1 - mu),   # Bernoulli
         "Gamma"             = mu^2,
         "inverse.gaussian"  = mu^3,
         { stop("Unsupported family for V(mu): ", family_name) }
  )
}

# Pearson dispersion: phi_hat = sum((y-mu)^2 / V(mu)) / (n - p)
# gaussian: phi = sigma^2; Gamma: shape = 1/phi; IG: lambda = 1/phi

#' Variance Function V(mu) for Common GLM Families
#'
#' Returns the mean–variance relationship \eqn{V(\mu)} used by 1-parameter
#' exponential-family GLMs. Implemented families:
#' \itemize{
#'   \item \code{"gaussian"}: \eqn{V(\mu)=1}
#'   \item \code{"poisson"}: \eqn{V(\mu)=\mu}
#'   \item \code{"binomial"} (Bernoulli): \eqn{V(\mu)=\mu(1-\mu)}
#'   \item \code{"Gamma"}: \eqn{V(\mu)=\mu^2}
#'   \item \code{"inverse.gaussian"}: \eqn{V(\mu)=\mu^3}
#' }
#'
#' @param mu Numeric vector of means.
#' @param family_name Character scalar; one of the supported family names
#'   (typically \code{family$family} from a GLM fit).
#'
#' @return Numeric vector \code{V(mu)} of the same length as \code{mu}.
#'
#' @examples
#' .V_mu(mu = c(1,2,3), family_name = "poisson")
#'
#' @keywords internal
#' @noRd
estimate_dispersion_pearson <- function(y, mu, X, family_name) {
  n <- length(y); p <- ncol(X)
  V <- .V_mu(mu, family_name)
  V[V <= .Machine$double.eps] <- .Machine$double.eps
  phi_hat <- sum((y - mu)^2 / V) / max(n - p, 1L)
  if (family_name %in% c("binomial","poisson")) return(1)  # canonical
  phi_hat
}

# Extract dispersion for a fitted margin from (y, X, eta, family)
# Returns: list(fam, link, disp, extra=list(size=...) for NB if present)

#' Extract/Estimate Dispersion (and NB Size) for a GLM Margin
#'
#' Estimates the dispersion parameter \eqn{\phi} for a fitted margin from
#' \code{(y, X, eta, family)}. For Negative Binomial, returns the \code{size}
#' (theta) instead. For Gamma, it attempts ML shape via \pkg{MASS}’s
#' \code{\link[MASS]{gamma.shape}}; if unavailable or it fails, it falls back
#' to the Pearson estimator.
#'
#' @param y Numeric response vector (\code{n}).
#' @param X Numeric design matrix (\code{n × p}); should include intercept column if modeled.
#' @param eta Numeric linear predictor \code{X*beta} (length \code{n}).
#' @param family A \code{\link[stats]{family}} object used for the margin.
#' @param fit_nb_size Optional numeric; if provided for Negative Binomial,
#'   it will be returned in \code{extra$size}.
#'
#' @details
#' \itemize{
#'   \item \strong{Gaussian}: returns \code{phi = sigma^2} (Pearson).
#'   \item \strong{Poisson/Binomial}: returns \code{phi = 1}.
#'   \item \strong{Gamma}: tries ML shape via \code{MASS::gamma.shape} if possible;
#'     otherwise Pearson. Final \code{phi = 1/shape}.
#'   \item \strong{Inverse-Gaussian}: Pearson; \code{lambda = 1/phi}.
#'   \item \strong{Negative Binomial}: returns \code{disp = 1} and
#'     \code{extra$size = theta}.
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{fam}}{Character family name.}
#'   \item{\code{link}}{Character link name.}
#'   \item{\code{disp}}{Estimated dispersion \eqn{\phi}.}
#'   \item{\code{extra}}{List with \code{size} when NB; otherwise empty.}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 200; X <- cbind(1, rnorm(n)); beta <- c(0, 0.5)
#' eta <- as.vector(X %*% beta); mu <- exp(eta)
#' y <- rgamma(n, shape = 2, scale = mu/2)
#' out <- extract_margin_dispersion(y, X, eta, family = Gamma("log"))
#' str(out)  # disp ≈ 1/shape
#'
#' @seealso \code{\link{estimate_dispersion_pearson}}, \code{\link[MASS]{gamma.shape}}
#' @export
extract_margin_dispersion <- function(y, X, eta, family, fit_nb_size = NULL) {
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

  # Negative Binomial handled via 'size' (theta)
  if (grepl("Negative Binomial", fam, ignore.case = TRUE)) {
    return(list(fam = "neg.binomial", link = link, disp = 1, extra = list(size = fit_nb_size)))
  }

  # --- Gamma: try ML shape via MASS::gamma.shape on a regular glm fit ---
  if (fam == "Gamma") {
    # Build a data.frame with expanded columns for X[,-1]
    if (requireNamespace("MASS", quietly = TRUE)) {
      # ensure strictly positive y for Gamma fit (defensive)
      y_pos <- pmax(y, .Machine$double.eps)

      if (ncol(X) > 1) {
        Wdf <- as.data.frame(X[, -1, drop = FALSE])
        colnames(Wdf) <- paste0("W", seq_len(ncol(Wdf)))
        df_tmp <- cbind.data.frame(y = y_pos, Wdf)
        form <- stats::as.formula(paste("y ~", paste(colnames(Wdf), collapse = " + ")))
      } else {
        df_tmp <- data.frame(y = y_pos)
        form <- y ~ 1
      }

      fit_tmp <- stats::glm(form, family = Gamma(link = "log"), data = df_tmp)
      sh <- tryCatch(MASS::gamma.shape(fit_tmp)$alpha, error = function(e) NA_real_)
      if (is.finite(sh) && sh > 0) {
        return(list(fam = fam, link = link, disp = 1/sh, extra = list()))
      }
      # fall through to Pearson if ML shape failed
    }

    # Pearson fallback if MASS not available or shape failed
    disp <- estimate_dispersion_pearson(y, mu, X, fam)
    return(list(fam = fam, link = link, disp = disp, extra = list()))
  }

  # --- Canonical 1-parameter EDFs ---
  if (fam %in% c("gaussian","poisson","binomial","inverse.gaussian")) {
    disp <- estimate_dispersion_pearson(y, mu, X, fam)
    return(list(fam = fam, link = link, disp = disp, extra = list()))
  }

  # Fallback (rare): variance of residuals
  disp <- stats::var(y - mu)
  list(fam = fam, link = link, disp = disp, extra = list())
}


# -------------------------------
# General CDF / PDF / Quantile
# -------------------------------


# MC fallback (rare)

#' Monte Carlo Fallback for GLM CDF (Internal)
#'
#' Crude Monte Carlo estimator of \eqn{F_Y(y|x)} used only as a last resort
#' when native or saddlepoint methods are not available/stable.
#' Implements simple random sampling for a few common families.
#'
#' @param y Numeric vector of evaluation points.
#' @param mu Numeric vector of conditional means (same length as \code{y} or length 1).
#' @param fam Character family name (e.g., \code{"gaussian"}, \code{"poisson"}, \code{"Gamma"}, \code{"binomial"}).
#' @param B Integer; number of Monte Carlo draws (default 20000).
#'
#' @return Numeric vector of approximate CDF values \code{P(Y ≤ y)}.
#'
#' @note This routine is intentionally conservative and not intended for
#' production inference; it’s only invoked when other methods fail.
#'
#' @examples
#' cdf_mc(y = c(0.5, 1, 2), mu = 1, fam = "poisson", B = 10000)
#'
#' @keywords internal
#' @noRd
cdf_mc <- function(y, mu, fam, B = 20000L) {
  y <- as.numeric(y); mu <- as.numeric(mu)
  if (length(mu)==1L && length(y)>1L) mu <- rep(mu, length(y))
  out <- numeric(length(y))
  for (i in seq_along(y)) {
    if (fam == "gaussian") {
      out[i] <- mean(stats::rnorm(B, mean = mu[i], sd = 1) <= y[i])
    } else if (fam == "poisson") {
      out[i] <- mean(stats::rpois(B, lambda = mu[i]) <= floor(y[i]))
    } else if (fam == "Gamma") {
      out[i] <- mean(stats::rgamma(B, shape = 1, scale = mu[i]) <= y[i]) # crude only
    } else if (fam == "binomial") {
      p <- max(min(mu[i],0.999),0.001)
      out[i] <- mean(stats::rbinom(B, size = 1, prob = p) <= y[i])
    } else {
      out[i] <- mean(stats::rnorm(B, mean = mu[i], sd = sqrt(abs(mu[i]))) <= y[i])
    }
  }
  out
}




# -------------------------------
# Core estimator (corrected)
# -------------------------------



# -------------------------------
# AB bootstrap wrapper
# -------------------------------


# =============================================================================
