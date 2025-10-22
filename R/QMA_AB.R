#' Adaptive Bootstrap Test for Quantile Mediation Effects
#'
#' Implements the adaptive bootstrap (AB) procedure from the paper,
#' combining the classical bootstrap and a local-limit correction
#' to handle the degenerate null (\eqn{\alpha=\beta=0}).
#'
#' @inheritParams estimate_quantile_mediation_AB
#' @param R Integer. Number of bootstrap replications.
#' @param lambda Numeric. Penalty scaling parameter; the test uses
#'   \eqn{\lambda_n = \lambda \sqrt{n}/\log n} (default 2).
#'
#' @return A numeric scalar p-value for testing the quantile indirect effect at \eqn{\tau}.
#'
#' @examples
#' dat <- gen_data_quantile_mediation()
#' W_cols <- 4:ncol(dat)
#' pval <- QMA_AB(
#'   data = dat,
#'   R = 100,
#'   tau = 0.5,
#'   x0 = 0, x1 = 1,
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
QMA_AB <- function(
    data,
    R = 500,
    lambda = 2,
    tau = 1/2,
    x0 = 0,
    x1 = 1,
    W_new = NULL,
    ...
) {
  n <- nrow(data)
  lambda_n <- lambda * sqrt(n) / log(n)

  boot_stat <- function(d, i) {
    est <- estimate_quantile_mediation_AB(
      data = d[i, , drop = FALSE],
      tau = tau, x0 = x0, x1 = x1, W_new = W_new, ...
    )
    unname(c(
      as.numeric(est$I_alpha_hat),
      as.numeric(est$I_beta_hat),
      est$NIQE_hat,
      est$alpha_hat,
      est$beta_hat,
      est$z_x0,
      est$z_x1,
      est$delta_Y_hat,
      est$eta_hat
    ))
  }

  b <- boot::boot(data = data, statistic = function(d,i) boot_stat(d,i), R = R)
  est <- boot_stat(data, seq_len(n))

  sd_n12_alpha <- sqrt(n) * stats::sd(b$t[, 4])
  sd_n12_beta  <- sqrt(n) * stats::sd(b$t[, 5])

  U <- apply(b$t, 1, function(x) {
    I_alpha_ast <- (abs(sqrt(n)*x[4]/sd_n12_alpha) <= lambda_n)
    I_beta_ast  <- (abs(sqrt(n)*x[5]/sd_n12_beta ) <= lambda_n)
    I_alpha     <- (abs(sqrt(n)*est[4]/sd_n12_alpha) <= lambda_n)
    I_beta      <- (abs(sqrt(n)*est[5]/sd_n12_beta ) <= lambda_n)
    indicator <- I_alpha_ast * I_alpha * I_beta_ast * I_beta

    (x[3] - est[3]) * (1 - indicator) +
      indicator * (1/n) *
      sqrt(n) * (x[4] - est[4]) *
      sqrt(n) * (x[5] - est[5]) *
      ((x[7] - x[6]) / x[8]) * x[9]
  })

  QNIE_hat <- est[3]
  mean(abs(U) > abs(QNIE_hat))
}
