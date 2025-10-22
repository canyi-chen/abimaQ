## ================================================================
## Parallel MC validation of p-value uniformity under the null
## ================================================================

# ---- libraries ----
library(boot)
library(mvtnorm)
library(future)
library(future.apply)
library(ggplot2)

# ---- single run to get a p-value under the null (alpha = 0) ----
one_pval <- function(n = 400, alpha = 0, beta = 0.8, gamma = 0.3,
                     tau = 0.5, Rboot = 400, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dat <- gen_data_quantile_mediation(n = n, alpha = alpha, beta = beta, gamma = gamma)
  W_cols <- 4:ncol(dat)
  W_new  <- matrix(0, 1, length(W_cols))
  QMA_AB(
    data = dat,
    R = Rboot,
    tau = tau,
    x0 = 0, x1 = 1,
    W_cols   = W_cols,
    W_X_cols = 1:2,
    W_M_cols = 1:length(W_cols),
    W_Y_cols = 1:length(W_cols),
    W_new    = W_new,
    fam_X = gaussian(),
    fam_M = gaussian(),
    fam_Y = Gamma(link = "log")
  )
}

# ---- simulation wrapper across many replicates and (optionally) taus ----
simulate_pvals <- function(Nrep = 1000, n = 400, tau = 0.5,
                           alpha = 0, beta = 0.8, gamma = 0.3,
                           Rboot = 400, parallel_workers = NULL) {
  if (is.null(parallel_workers)) {
    parallel_workers <- max(1, parallel::detectCores() - 1)
  }
  plan(multisession, workers = parallel_workers)

  # future.apply ensures independent RNG streams if future.seed=TRUE
  pvals <- future_sapply(
    X = seq_len(Nrep),
    FUN = function(i) one_pval(n = n, alpha = alpha, beta = beta, gamma = gamma,
                               tau = tau, Rboot = Rboot),
    future.seed = TRUE
  )

  # basic diagnostics
  rates <- c(
    "alpha_1%"  = mean(pvals < 0.01),
    "alpha_5%"  = mean(pvals < 0.05),
    "alpha_10%" = mean(pvals < 0.10)
  )
  ks <- suppressWarnings(ks.test(pvals, "punif"))

  list(
    pvals = pvals,
    rates = rates,
    ks    = ks,
    meta  = list(Nrep = Nrep, n = n, tau = tau,
                 alpha = alpha, beta = beta, gamma = gamma, Rboot = Rboot)
  )
}

# ---- plotting helpers ----
plot_pval_diagnostics <- function(pvals, title_suffix = "") {
  df <- data.frame(p = pvals)
  N  <- length(pvals)

  # Histogram with uniform overlay
  g1 <- ggplot(df, aes(x = p)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, boundary = 0, closed = "left") +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(title = paste0("P-value Histogram", title_suffix),
         x = "p", y = "density") +
    ggpubr::theme_pubr()

  # ECDF vs diagonal
  g2 <- ggplot(df, aes(x = p)) +
    stat_ecdf(geom = "step") +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_equal() +
    labs(title = paste0("ECDF of p vs Uniform(0,1)", title_suffix),
         x = "p", y = "F_n(p)") +
    ggpubr::theme_pubr()

  # QQ vs Uniform(0,1)
  p_sorted <- sort(pvals)
  u_theory <- (seq_len(N) - 0.5) / N
  dfqq <- data.frame(theory = u_theory, sample = p_sorted)
  g3 <- ggplot(dfqq, aes(x = theory, y = sample)) +
    geom_point(size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_equal() +
    labs(title = paste0("QQ plot of p values", title_suffix),
         x = "Theoretical quantiles", y = "Empirical quantiles") +
    ggpubr::theme_pubr()

  list(hist = g1, ecdf = g2, qq = g3)
}

if(FALSE) {
  # ================================================================
  # RUN: (adjust Nrep to 200–1000 depending on runtime budget)
  # ================================================================
  set.seed(2025)
  res_null <- simulate_pvals(
    Nrep = 500,       # bump to 1000+ for a thorough check
    n    = 400,
    tau  = 0.5,
    alpha = 0,        # NULL via alpha=0 (so alpha*beta=0)
    beta  = 1,
    gamma = 1,
    Rboot = 300       # increase for higher-precision per-rep p-values
  )

  print(res_null$rates)       # rejection rates at 1%, 5%, 10%
  print(res_null$ks)          # KS test vs Uniform(0,1)

  plots <- plot_pval_diagnostics(res_null$pvals)
  print(plots$hist)
  print(plots$ecdf)
  print(plots$qq)

  # ---- Optional: other corner of the null (beta = 0) ----
  # res_null_beta0 <- simulate_pvals(Nrep = 500, n = 400, tau = 0.5,
  #                                  alpha = 0.6, beta = 0, gamma = 0.3, Rboot = 300)
  # print(res_null_beta0$rates)
  # print(res_null_beta0$ks)
  # plots2 <- plot_pval_diagnostics(res_null_beta0$pvals, title_suffix = " (beta=0 null)")
  # print(plots2$hist); print(plots2$ecdf); print(plots2$qq)

  # ---- Clean up parallel workers when done ----
  plan(sequential)
}
