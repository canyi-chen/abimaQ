#' Estimate Quantile Mediation Effects (Gaussian Copula + GLM Margins)
#'
#' Fits a Gaussian-copula structural equation model with GLM margins for
#' exposure (\eqn{X}), mediator (\eqn{M}), and outcome (\eqn{Y}),
#' and computes plug-in estimators for quantile natural effects (Theorem 1).
#'
#' @param data Data frame with columns \code{X}, \code{M}, \code{Y}, then covariates \code{W}.
#' @param tau Numeric, quantile level \eqn{\tau \in (0,1)}.
#' @param x0,x1 Numeric, counterfactual exposure levels.
#' @param W_cols Integer indices of \code{W} columns in \code{data}.
#' @param W_X_cols,W_M_cols,W_Y_cols Integer indices (within \code{W_cols})
#'   used in the \eqn{X|W}, \eqn{M|W}, and \eqn{Y|W} margins, respectively.
#' @param W_new 1 \eqn{\times} p matrix of covariate values for counterfactual evaluation.
#' @param fam_X,fam_M,fam_Y \code{\link[stats]{family}} objects for the three margins.
#' @param disp_X,disp_M,disp_Y Numeric dispersion parameters (estimated internally if omitted).
#' @param nb_size Optional size (theta) for Negative-Binomial margins.
#' @param pit_mode One of \code{"midp"}, \code{"jitter"}, \code{"closed"} for discrete PIT handling.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{NIQE_hat}}{Estimated quantile indirect effect at \eqn{\tau}.}
#'   \item{\code{alpha_hat}, \code{beta_hat}, \code{gamma_hat}}{Estimated copula parameters.}
#'   \item{\code{delta_Y_hat}}{Scale factor for the latent \eqn{Y}.}
#'   \item{\code{z_x0}, \code{z_x1}}{Probit PITs for \eqn{X} at \eqn{x_0}, \eqn{x_1}.}
#'   \item{\code{I_alpha_hat}, \code{I_beta_hat}}{Pretest indicators.}
#' }
#'
#' @examples
#' dat <- gen_data_quantile_mediation()
#' W_cols <- 4:ncol(dat)
#' est <- estimate_quantile_mediation_AB(
#'   data = dat,
#'   tau = 0.5,
#'   x0 = 0,
#'   x1 = 1,
#'   W_cols = W_cols,
#'   W_X_cols = 1:2,
#'   W_M_cols = 1:length(W_cols),
#'   W_Y_cols = 1:length(W_cols),
#'   W_new = matrix(0, 1, length(W_cols)),
#'   fam_X = gaussian(),
#'   fam_M = gaussian(),
#'   fam_Y = Gamma(link = "log")
#' )
#' @export
estimate_quantile_mediation_AB <- function(
    data,
    tau = 1/2,
    x0 = 0,
    x1 = 1,
    W_new = NULL,
    W_cols = NULL,
    W_X_cols = NULL,
    W_M_cols = NULL,
    W_Y_cols = NULL,
    fam_X = gaussian(),
    fam_M = gaussian(),
    fam_Y = Gamma(link="log"),
    disp_X = 1,
    disp_M = 1,
    disp_Y = 1,
    nb_size = NULL,
    pit_mode = c("midp","jitter","closed")
) {
  pit_mode <- match.arg(pit_mode)
  stopifnot(tau > 0 && tau < 1)

  data <- as.data.frame(data)
  n <- nrow(data); if (n < 5L) stop("n too small.")
  X <- as.numeric(data[[1]])
  M <- as.numeric(data[[2]])
  Y <- as.numeric(data[[3]])

  if (is.null(W_cols)) {
    if (ncol(data) < 4L) stop("No W columns detected; set W_cols or include W in data.")
    W_cols <- 4:ncol(data)
  }
  W_all <- as.matrix(data[, W_cols, drop = FALSE])
  pW <- ncol(W_all)
  if (is.null(W_new)) W_new <- matrix(0, 1, pW)
  W_new <- as.matrix(W_new)
  if (ncol(W_new) != pW) stop("W_new must have ncol = length(W_cols).")

  if (is.null(W_X_cols)) W_X_cols <- seq_len(pW)
  if (is.null(W_M_cols)) W_M_cols <- seq_len(pW)
  if (is.null(W_Y_cols)) W_Y_cols <- seq_len(pW)

  W_X <- W_all[, W_X_cols, drop = FALSE]
  W_M <- W_all[, W_M_cols, drop = FALSE]
  W_Y <- W_all[, W_Y_cols, drop = FALSE]

  WX_new <- W_new[, W_X_cols, drop = FALSE]
  WY_new <- W_new[, W_Y_cols, drop = FALSE]

  ones <- rep(1, n)
  XWX <- cbind(ones, W_X)
  XWM <- cbind(ones, W_M)
  XWY <- cbind(ones, W_Y)

  XWX_new <- cbind(1, WX_new)
  XWY_new <- cbind(1, WY_new)

  use_fastglm <- function() requireNamespace("fastglm", quietly = TRUE)
  fit_gaussian_sigma <- function(X, y) {
    QR <- .lm.fit(X, y)
    mu <- drop(X %*% QR$coefficients)
    p  <- ncol(X)
    rss <- sum((y - mu)^2)
    sig2 <- rss / pmax(n - p, 1)
    list(coef = QR$coefficients, mu = mu, sigma2 = sig2)
  }
  fit_glm <- function(X, y, family) {
    famname <- family$family
    if (grepl("Negative Binomial", famname, ignore.case = TRUE)) {
      if (!requireNamespace("MASS", quietly = TRUE))
        stop("MASS not available for Negative Binomial fit.")
      fit <- MASS::glm.nb(y ~ X - 1)
      return(list(coef = stats::coef(fit), family = "neg.binomial", link = fit$family$link,
                  nb_size = fit$theta, fit = fit))
    }
    if (use_fastglm() && !(family$family %in% c("binomial"))) {
      co <- fastglm::fastglm(X, y, family = family)$coefficients
    } else {
      gf <- stats::glm.fit(x = X, y = y, family = family)
      co <- gf$coefficients
    }
    list(coef = co, family = family$family, link = family$link)
  }

  # X|W_X
  if (identical(fam_X$family, gaussian()$family)) {
    fx <- fit_gaussian_sigma(XWX, X)
    eta_X <- drop(XWX %*% fx$coef)
    mu_X  <- fx$mu
    disp_X <- fx$sigma2
    coef_X <- fx$coef
  } else {
    fx <- fit_glm(XWX, X, fam_X)
    eta_X <- drop(XWX %*% fx$coef)
    mu_X  <- switch(fx$link,
                    "identity"=eta_X,"log"=exp(eta_X),"inverse"=1/eta_X,
                    "logit"=1/(1+exp(-eta_X)),"probit"=stats::pnorm(eta_X),
                    "cloglog"=1-exp(-exp(eta_X)), exp(eta_X))
    coef_X <- fx$coef
  }
  disp_info_X <- extract_margin_dispersion(y = X, X = XWX, eta = eta_X, family = fam_X)
  disp_X <- disp_info_X$disp

  # M|W_M
  if (identical(fam_M$family, gaussian()$family)) {
    fm <- fit_gaussian_sigma(XWM, M)
    eta_M <- drop(XWM %*% fm$coef)
    mu_M  <- fm$mu
    disp_M <- fm$sigma2
    coef_M <- fm$coef
  } else {
    fm <- fit_glm(XWM, M, fam_M)
    eta_M <- drop(XWM %*% fm$coef)
    mu_M  <- switch(fm$link,
                    "identity"=eta_M,"log"=exp(eta_M),"inverse"=1/eta_M,
                    "logit"=1/(1+exp(-eta_M)),"probit"=stats::pnorm(eta_M),
                    "cloglog"=1-exp(-exp(eta_M)), exp(eta_M))
    coef_M <- fm$coef
  }
  disp_info_M <- extract_margin_dispersion(y = M, X = XWM, eta = eta_M, family = fam_M)
  disp_M <- disp_info_M$disp

  # Y|W_Y
  fy <- fit_glm(XWY, Y, fam_Y)
  eta_Y <- drop(XWY %*% fy$coef)
  mu_Y  <- switch(fy$link,
                  "identity"=eta_Y,"log"=exp(eta_Y),"inverse"=1/eta_Y,
                  "logit"=1/(1+exp(-eta_Y)),"probit"=stats::pnorm(eta_Y),
                  "cloglog"=1-exp(-exp(eta_Y)), exp(eta_Y))
  coef_Y <- fy$coef
  if (grepl("Negative Binomial", fam_Y$family, ignore.case = TRUE) && is.null(nb_size)) {
    if (!is.null(fy$nb_size)) nb_size <- fy$nb_size
  }
  disp_info_Y <- extract_margin_dispersion(y = Y, X = XWY, eta = eta_Y, family = fam_Y,
                                           fit_nb_size = nb_size)
  disp_Y <- disp_info_Y$disp
  if (is.null(nb_size) && length(disp_info_Y$extra$size)) nb_size <- disp_info_Y$extra$size

  # PIT -> Z (stable; no row dropping)
  eps <- 1e-12; clip01 <- function(p) pmin(pmax(p, eps), 1 - eps)
  Z_X <- stats::qnorm( clip01(cdf_glm_general(X, eta = eta_X, family = fam_X,
                                              dispersion = disp_X, method = "auto",
                                              discrete_pit = pit_mode, nb_size = nb_size)) )
  Z_M <- stats::qnorm( clip01(cdf_glm_general(M, eta = eta_M, family = fam_M,
                                              dispersion = disp_M, method = "auto",
                                              discrete_pit = pit_mode, nb_size = nb_size)) )
  Z_Y <- stats::qnorm( clip01(cdf_glm_general(Y, eta = eta_Y, family = fam_Y,
                                              dispersion = disp_Y, method = "auto",
                                              discrete_pit = pit_mode, nb_size = nb_size)) )
  Z <- cbind(Z_X, Z_M, Z_Y)
  n_eff <- nrow(Z)

  # Copula likelihood for Z ~ N(0, Sigma(theta))
  nmloglik_mvn0 <- function(Z, Sigma) {
    R <- chol(Sigma); logdet <- 2*sum(log(diag(R)))
    Yt <- forwardsolve(t(R), t(Z)); quad <- colSums(Yt^2)
    0.5*(logdet + mean(quad) + ncol(Z)*log(2*pi))
  }
  loss <- function(theta) {
    a <- theta[1]; b <- theta[2]; g <- theta[3]
    e <- g + a*b
    dM <- sqrt(a^2 + 1); dY <- sqrt(e^2 + b^2 + 1)
    s12 <- a/dM; s13 <- e/dY; s23 <- (a*e + b)/(dM*dY)
    Sigma <- matrix(c(1,s12,s13,s12,1,s23,s13,s23,1),3,3,byrow=TRUE)
    nmloglik_mvn0(Z, Sigma)
  }
  opt <- optim(c(0,0,0), fn = loss, method = "BFGS",
               control = list(reltol = 1e-10, maxit = 1000))
  alpha_hat <- opt$par[1]; beta_hat <- opt$par[2]; gamma_hat <- opt$par[3]
  delta_Y_hat <- sqrt((gamma_hat + alpha_hat*beta_hat)^2 + beta_hat^2 + 1)

  # Counterfactuals at W_new
  etaX_new <- drop(cbind(1, WX_new) %*% coef_X)
  z_x0 <- stats::qnorm( clip01(cdf_glm_general(x0, eta=etaX_new, family=fam_X,
                                               dispersion=disp_X, method="auto",
                                               discrete_pit=pit_mode, nb_size=nb_size)) )
  z_x1 <- stats::qnorm( clip01(cdf_glm_general(x1, eta=etaX_new, family=fam_X,
                                               dispersion=disp_X, method="auto",
                                               discrete_pit=pit_mode, nb_size=nb_size)) )

  etaY_new <- drop(cbind(1, WY_new) %*% coef_Y)

  transform_p <- function(z0, z1) {
    stats::pnorm( (gamma_hat*z0 + alpha_hat*beta_hat*z1 + stats::qnorm(tau)*sqrt(1+beta_hat^2)) / delta_Y_hat )
  }
  p1 <- transform_p(z_x1, z_x1)
  p0 <- transform_p(z_x1, z_x0)

  Fy_inv <- function(p) q_glm_general(p, eta = etaY_new, family = fam_Y,
                                      dispersion = disp_Y, nb_size = nb_size)
  fY <- function(y) pdf_glm_general(y, eta = etaY_new, family = fam_Y,
                                    dispersion = disp_Y, method = "auto", nb_size = nb_size)

  NIQE_hat <- Fy_inv(p1) - Fy_inv(p0)

  I_alpha_hat <- (abs(sqrt(n_eff)*alpha_hat) <= 2*sqrt(n_eff)/log(n_eff))
  I_beta_hat  <- (abs(sqrt(n_eff)*beta_hat ) <= 2*sqrt(n_eff)/log(n_eff))

  w_hat  <- (gamma_hat*z_x1 + stats::qnorm(tau)*sqrt(1 + beta_hat^2)) / delta_Y_hat
  qy     <- Fy_inv(stats::pnorm(w_hat))
  eta_hat <- stats::dnorm(w_hat) / fY(qy)

  list(
    I_alpha_hat = I_alpha_hat,
    I_beta_hat  = I_beta_hat,
    NIQE_hat    = as.numeric(NIQE_hat),
    alpha_hat   = alpha_hat,
    beta_hat    = beta_hat,
    z_x0        = z_x0,
    z_x1        = z_x1,
    delta_Y_hat = delta_Y_hat,
    eta_hat     = eta_hat
  )
}
