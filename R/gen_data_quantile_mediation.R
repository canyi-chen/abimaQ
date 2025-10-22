#' @title Generate Synthetic Data for Quantile Mediation Analysis
#' @description
#' Simulates data under a Gaussian copula structural equation model with
#' Gaussian margins for exposure (X) and mediator (M), and Gamma (Exponential) margins for outcome (Y).
#' The design matches the data generating process in the paper for validating the AB test under the null.
#'
#' @param n Sample size (default `200`).
#' @param alpha Path coefficient from X to M on the latent scale.
#' @param beta Path coefficient from M to Y on the latent scale.
#' @param gamma Path coefficient from X to Y on the latent scale.
#' @param beta_X Coefficients for W in the X model.
#' @param beta_M Coefficients for W in the M model.
#' @param beta_Y Coefficients for W in the Y model.
#' @param sigmas Standard deviations for X and M margins.
#' @param rho_W Compound symmetry correlation for W.
#'
#' @return A data frame with columns \code{X}, \code{M}, \code{Y} and \code{W}.
#' @examples
#' set.seed(1)
#' dat <- gen_data_quantile_mediation(n = 500, alpha = 0, beta = 1, gamma = 1)
#' head(dat)
#' @export
gen_data_quantile_mediation <- function(n = 2e2,
                                        alpha = 0,
                                        beta = 1,
                                        gamma = 1,
                                        beta_X = c(0.5, 0.2, 0.2),
                                        beta_M = c(0.8, 0.3, 0.3, 0.4),
                                        beta_Y = c(-0.2, 0.4, -0.2, 0.7),
                                        sigmas = rep(0.3, 2),
                                        rho_W = 0.2
                                        # sigmas = rep(1, 3) # too noisy
) {

  p1 <- length(beta_X)
  p <- length(beta_M)
  p2 <- p - p1


  # generate W1 and W from a multivariate normal
  # with mean zero and compound symmetry correlation
  Sigma_CS <- matrix(rho_W, p - 1, p - 1)
  diag(Sigma_CS) <- 1
  W <- mvtnorm::rmvnorm(n, rep(0, p - 1), sigma = Sigma_CS)
  W <- cbind(1, W)
  W1 <- W[, 1:p1]


  Z_X <- rnorm(n)
  Z_M <- rnorm(n, mean = alpha * Z_X)
  Z_Y <- rnorm(n, mean = gamma * Z_X + beta * Z_M)


  delta_M <- sqrt(alpha^2 + 1)
  delta_Y <- sqrt((gamma + alpha * beta)^2 + beta^2 + 1)
  Z_M_ast <- Z_M / delta_M
  Z_Y_ast <- Z_Y / delta_Y

  X <- qnorm(pnorm(Z_X), mean = W1 %*% beta_X, sd = sigmas[1])
  M <- qnorm(pnorm(Z_M_ast), mean = W %*% beta_M, sd = sigmas[2])
  Y <- qgamma(pnorm(Z_Y_ast), shape = 1, scale = exp(W %*% beta_Y))
  return(data.frame(X = X,
                    M = M,
                    Y = Y,
                    W = W[, -1]))
}
# --- your generator (as given) ---
gen_data_quantile_mediation <- function(n = 2e2,
                                        alpha = 0,
                                        beta = 1,
                                        gamma = 1,
                                        beta_X = c(0.5, 0.2, 0.2),
                                        beta_M = c(0.8, 0.3, 0.3, 0.4),
                                        beta_Y = c(-0.2, 0.4, -0.2, 0.7),
                                        sigmas = rep(0.3, 2),
                                        rho_W = 0.2) {
  p1 <- length(beta_X)
  p  <- length(beta_M)
  p2 <- p - p1

  # W ~ MVN with compound symmetry (needs mvtnorm)
  Sigma_CS <- matrix(rho_W, p - 1, p - 1); diag(Sigma_CS) <- 1
  W <- mvtnorm::rmvnorm(n, rep(0, p - 1), sigma = Sigma_CS)
  W <- cbind(1, W) # add intercept
  W1 <- W[, 1:p1]

  # latent normals with target correlations
  Z_X <- rnorm(n)
  Z_M <- rnorm(n, mean = alpha * Z_X)
  Z_Y <- rnorm(n, mean = gamma * Z_X + beta * Z_M)

  delta_M <- sqrt(alpha^2 + 1)
  delta_Y <- sqrt((gamma + alpha * beta)^2 + beta^2 + 1)
  Z_M_ast <- Z_M / delta_M
  Z_Y_ast <- Z_Y / delta_Y

  # observed margins
  X <- qnorm(pnorm(Z_X),     mean = W1 %*% beta_X, sd = sigmas[1])
  M <- qnorm(pnorm(Z_M_ast), mean = W %*% beta_M,  sd = sigmas[2])
  Y <- qgamma(pnorm(Z_Y_ast), shape = 1, scale = exp(W %*% beta_Y))  # Exp with mean exp(W%*%beta_Y)

  data.frame(X = X, M = M, Y = Y, W = W[, -1]) # drop intercept in returned W columns
}

# --- demo: generate data, run estimator & AB test (uses functions I sent) ---
demo_qma <- function(n = 400, alpha = 0.6, beta = 0.8, gamma = 0.3,
                     tau = 0.5, R = 200, seed = 1) {
  set.seed(seed)
  dat <- gen_data_quantile_mediation(n = n, alpha = alpha, beta = beta, gamma = gamma)

  # figure out W columns automatically (X,M,Y are first 3)
  W_cols <- 4:ncol(dat)

  # one new W point to evaluate counterfactuals; 0-vector (same length as W_cols)
  W_new <- matrix(0, 1, length(W_cols))

  est <- estimate_quantile_mediation_AB(
    data = dat,
    tau = tau,
    x0 = 0, x1 = 1,
    W_cols = W_cols,
    W_X_cols = 1:2,          # use first two W's in X|W (tweak as you like)
    W_M_cols = 1:length(W_cols),  # use all W's in M|W
    W_Y_cols = 1:length(W_cols),  # use all W's in Y|W
    W_new = W_new,
    fam_X = gaussian(),
    fam_M = gaussian(),
    fam_Y = Gamma(link = "log")   # matches the generator (shape=1 -> φ=1)
  )

  pval <- QMA_AB(
    data = dat,
    R = R,
    tau = tau,
    x0 = 0, x1 = 1,
    W_cols = W_cols,
    W_X_cols = 1:2,
    W_M_cols = 1:length(W_cols),
    W_Y_cols = 1:length(W_cols),
    W_new = W_new,
    fam_X = gaussian(),
    fam_M = gaussian(),
    fam_Y = Gamma(link = "log")
  )

  list(est = est, pval = pval)
}

# Example run:
# out <- demo_qma()
# str(out)
# out$pval
# out$est$NIQE_hat
