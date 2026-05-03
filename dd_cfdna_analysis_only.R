library(rjags)
library(loo)
library(pbivnorm)
library(readxl)
library(coda)
library(mada)
library(MASS)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(tidyr)

rho_values  <- c(0.00, 0.30, 0.50, 0.70)
primary_rho <- 0.50

stopifnot(primary_rho %in% rho_values)

rho_repeat_default <- 0.50

mcmc_main <- list(adapt = 5000, burnin = 25000, iter = 100000,
                  thin = 10, chains = 4)
mcmc_loo  <- list(adapt = 2000, burnin = 10000, iter = 50000,
                  thin = 10, chains = 2)
mcmc_quad <- list(adapt = 3000, burnin = 15000, iter = 50000,
                  thin = 10, chains = 2)

SEED_MASTER <- 42L

DNA <- read_excel("DNA.xlsx")

df_raw <- DNA %>%
  rename(
    study     = Study,
    threshold = Theshold,
    TP = TP, FP = FP, FN = FN, TN = TN
  ) %>%
  mutate(
    n_pos = TP + FN,
    n_neg = FP + TN,
    sens  = TP / n_pos,
    spec  = TN / n_neg,
    fpr   = 1 - spec,

    threshold_c = threshold - mean(threshold)
  )

disambiguate_studies <- function(study_vec) {
  monica_idx <- which(study_vec == "Monica et al.")
  if (length(monica_idx) == 2) {
    study_vec[monica_idx[1]] <- "Monica et al. (a)"
    study_vec[monica_idx[2]] <- "Monica et al. (b)"
  }
  study_vec
}
df_raw$study <- disambiguate_studies(df_raw$study)

cat("\nStudy labels after disambiguation:\n")
for (i in seq_len(nrow(df_raw))) {
  cat(sprintf("  %2d. %s\n", i, df_raw$study[i]))
}
n_unique <- length(unique(df_raw$study))
cat(sprintf("Unique study names: %d / %d\n", n_unique, nrow(df_raw)))
if (n_unique != nrow(df_raw)) {
  warning("Duplicate study labels remain after disambiguation: ",
          paste(df_raw$study[duplicated(df_raw$study)], collapse = ", "),
          call. = FALSE)
}

pat_col <- intersect(c("No", "Patients", "N_patients"), names(DNA))
if (length(pat_col) == 0) {
  stop("No patient-count column found in DNA.xlsx. ",
       "The cluster-correction DEff requires a per-study patient count. ",
       "Add a column named 'No', 'Patients', or 'N_patients'.")
}
df_raw$n_patients <- DNA[[pat_col[1]]]
df_raw <- df_raw %>%
  mutate(
    n_tests = n_pos + n_neg,
    m_bar   = n_tests / n_patients
  )

N <- nrow(df_raw)
mean_thresh <- mean(df_raw$threshold)

cat(sprintf("\nStudies: %d | Threshold range: %.3f-%.3f (mean = %.3f)\n",
            N, min(df_raw$threshold), max(df_raw$threshold), mean_thresh))
cat(sprintf("Patient-count column: %s\n", pat_col[1]))

apply_deff <- function(df, rho) {
  preserve_zero <- function(x, deff) {
    out <- ifelse(x == 0, 0, pmax(1, round(x / deff)))
    as.integer(out)
  }
  df %>%
    mutate(
      deff      = 1 + (m_bar - 1) * rho,
      TP_eff    = preserve_zero(TP, deff),
      FP_eff    = preserve_zero(FP, deff),
      FN_eff    = preserve_zero(FN, deff),
      TN_eff    = preserve_zero(TN, deff),
      n_pos_eff = TP_eff + FN_eff,
      n_neg_eff = FP_eff + TN_eff
    )
}

cat("\n=== DEff Sensitivity Table (cluster-corrected: 1 + (m_bar - 1)*rho) ===\n")
deff_summary <- lapply(rho_values, function(rho) {
  d <- apply_deff(df_raw, rho)
  data.frame(
    rho            = rho,
    mean_deff      = mean(d$deff),
    median_deff    = median(d$deff),
    min_deff       = min(d$deff),
    max_deff       = max(d$deff),
    total_n_orig   = sum(d$n_pos + d$n_neg),
    total_n_eff    = sum(d$n_pos_eff + d$n_neg_eff),
    pct_retained   = 100 * sum(d$n_pos_eff + d$n_neg_eff) /
                     sum(d$n_pos + d$n_neg)
  )
}) %>% do.call(rbind, .)
print(deff_summary, row.names = FALSE, digits = 3)

cat("\nPer-study DEff (each column is a rho value):\n")
deff_per_study <- df_raw %>%
  transmute(study, n_tests, n_patients, m_bar = round(m_bar, 2)) %>%
  as.data.frame()
for (rho in rho_values) {
  deff_per_study[[sprintf("DEff_rho_%.2f", rho)]] <-
    round(1 + (df_raw$m_bar - 1) * rho, 2)
}
print(deff_per_study, row.names = FALSE)
cat(sprintf("\n* Primary analysis: rho = %.2f\n", primary_rho))
cat("  rho = 1.00 (not in sweep) reproduces the OLD DEff = n_tests/n_patients formula.\n")

hsroc_model_string <- "
model {
  for (i in 1:N) {
    lambda[i] <- gamma_0 + gamma_1 * threshold_c[i] + u[i]
    u[i]       ~ dnorm(0, tau_u)

    theta[i]   ~ dnorm(mu_theta, tau_theta)

    alpha[i] <- theta[i]  * exp(-beta / 2)
    psi[i]   <- lambda[i] * exp( beta / 2)

    logit_sens[i] <- (alpha[i] + psi[i]) / 2
    logit_fpr[i]  <- (-alpha[i] + psi[i]) / 2

    p_sens[i] <- exp(logit_sens[i]) / (1 + exp(logit_sens[i]))
    p_fpr[i]  <- exp(logit_fpr[i])  / (1 + exp(logit_fpr[i]))

    TP[i] ~ dbin(p_sens[i], n_pos[i])
    FP[i] ~ dbin(p_fpr[i],  n_neg[i])

    # Per-study log-likelihood for PSIS-LOO (joint TP + FP contribution).
    # Uses the study-specific p_sens[i], p_fpr[i] from the full-data fit -
    # this is the standard log_lik that loo::loo() expects.
    log_lik[i] <- logdensity.bin(TP[i], p_sens[i], n_pos[i]) +
                  logdensity.bin(FP[i], p_fpr[i], n_neg[i])
  }

  mu_theta  ~ dnorm(0, 1)
  gamma_0   ~ dnorm(0, 0.25)
  gamma_1   ~ dnorm(0, 0.01)
  beta      ~ dnorm(0, 4)

  tau_theta ~ dgamma(2, 0.5)
  tau_u     ~ dgamma(2, 0.5)

  sigma_theta <- 1 / sqrt(tau_theta)
  sigma_u     <- 1 / sqrt(tau_u)

  # Posterior predictive SROC over a fine dd-cfDNA % grid (population mean)
  for (g in 1:G) {
    lambda_grid[g] <- gamma_0 + gamma_1 * t_grid_c[g]
    alpha_grid[g]  <- mu_theta       * exp(-beta / 2)
    psi_grid[g]    <- lambda_grid[g] * exp( beta / 2)

    logit_s_grid[g] <- (alpha_grid[g] + psi_grid[g]) / 2
    logit_f_grid[g] <- (-alpha_grid[g] + psi_grid[g]) / 2

    sens_grid[g] <- exp(logit_s_grid[g]) / (1 + exp(logit_s_grid[g]))
    fpr_grid[g]  <- exp(logit_f_grid[g]) / (1 + exp(logit_f_grid[g]))
    spec_grid[g] <- 1 - fpr_grid[g]
  }
}
"

hsroc_quad_string <- "
model {
  for (i in 1:N) {
    lambda[i] <- gamma_0 + gamma_1 * threshold_c[i] + gamma_2 * threshold_c2[i] + u[i]
    u[i]       ~ dnorm(0, tau_u)
    theta[i]   ~ dnorm(mu_theta, tau_theta)

    alpha[i] <- theta[i]  * exp(-beta / 2)
    psi[i]   <- lambda[i] * exp( beta / 2)

    logit_sens[i] <- (alpha[i] + psi[i]) / 2
    logit_fpr[i]  <- (-alpha[i] + psi[i]) / 2

    p_sens[i] <- exp(logit_sens[i]) / (1 + exp(logit_sens[i]))
    p_fpr[i]  <- exp(logit_fpr[i])  / (1 + exp(logit_fpr[i]))

    TP[i] ~ dbin(p_sens[i], n_pos[i])
    FP[i] ~ dbin(p_fpr[i],  n_neg[i])
  }

  mu_theta  ~ dnorm(0, 1)
  gamma_0   ~ dnorm(0, 0.25)
  gamma_1   ~ dnorm(0, 0.01)
  gamma_2   ~ dnorm(0, 0.01)
  beta      ~ dnorm(0, 4)
  tau_theta ~ dgamma(2, 0.5)
  tau_u     ~ dgamma(2, 0.5)
  sigma_theta <- 1 / sqrt(tau_theta)
  sigma_u     <- 1 / sqrt(tau_u)
}
"

min_t      <- max(min(df_raw$threshold) * 0.5, 0.01)
max_t      <- max(df_raw$threshold) * 1.25
t_grid_raw <- seq(min_t, max_t, length.out = 200)
t_grid_c   <- t_grid_raw - mean_thresh
G          <- length(t_grid_raw)

cat(sprintf("\nThreshold grid: %.3f%%-%.3f%% (max observed = %.3f%%, headroom = %.3f%%)\n",
            min_t * 100, max_t * 100, max(df_raw$threshold) * 100, max_t * 100))

make_inits <- function(n_chains, base_seed) {
  lapply(seq_len(n_chains), function(i) {
    list(
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = as.integer(base_seed + i - 1L)
    )
  })
}

exact_loo_marginal <- function(i, df_eff, hsroc_model_string,
                                t_grid_c, G, mcmc_loo,
                                seed_master, rho_index,
                                n_inner = 1000L) {
  jags_data_loo <- list(
    N           = nrow(df_eff) - 1L,
    TP          = df_eff$TP_eff[-i],
    FP          = df_eff$FP_eff[-i],
    n_pos       = df_eff$n_pos_eff[-i],
    n_neg       = df_eff$n_neg_eff[-i],
    threshold_c = df_eff$threshold_c[-i],
    G           = G,
    t_grid_c    = t_grid_c
  )
  inits_loo <- make_inits(mcmc_loo$chains,
                           base_seed = seed_master + rho_index * 100000L + i)
  m <- jags.model(textConnection(hsroc_model_string),
                   data = jags_data_loo,
                   inits = inits_loo,
                   n.chains = mcmc_loo$chains,
                   n.adapt = mcmc_loo$adapt, quiet = TRUE)
  update(m, n.iter = mcmc_loo$burnin, progress.bar = "none")
  loo_post <- coda.samples(m,
    variable.names = c("mu_theta", "gamma_0", "gamma_1", "beta",
                       "sigma_theta", "sigma_u"),
    n.iter = mcmc_loo$iter, thin = mcmc_loo$thin,
    progress.bar = "none")
  ld <- do.call(rbind, loo_post)
  D  <- nrow(ld)

  tc_i   <- df_eff$threshold_c[i]
  tp_i   <- df_eff$TP_eff[i];     fp_i   <- df_eff$FP_eff[i]
  npos_i <- df_eff$n_pos_eff[i];  nneg_i <- df_eff$n_neg_eff[i]

  set.seed(seed_master + rho_index * 1000000L + i * 31L)
  log_q <- numeric(D)
  for (d in seq_len(D)) {
    mu  <- ld[d, "mu_theta"]; g0 <- ld[d, "gamma_0"]
    g1  <- ld[d, "gamma_1"];  bd <- ld[d, "beta"]
    sth <- ld[d, "sigma_theta"]; su <- ld[d, "sigma_u"]

    th_m <- rnorm(n_inner, mu, sth)
    u_m  <- rnorm(n_inner, 0, su)
    la_m <- g0 + g1 * tc_i + u_m
    al_m <- th_m * exp(-bd / 2)
    ps_m <- la_m * exp( bd / 2)

    p_sens_m <- plogis((al_m + ps_m) / 2)
    p_fpr_m  <- plogis((-al_m + ps_m) / 2)
    p_sens_m <- pmin(pmax(p_sens_m, 1e-10), 1 - 1e-10)
    p_fpr_m  <- pmin(pmax(p_fpr_m,  1e-10), 1 - 1e-10)

    log_lik_inner <- dbinom(tp_i, npos_i, p_sens_m, log = TRUE) +
                     dbinom(fp_i, nneg_i, p_fpr_m,  log = TRUE)
    mll_inner <- max(log_lik_inner)
    log_q[d]  <- mll_inner + log(mean(exp(log_lik_inner - mll_inner)))
  }
  mll_outer <- max(log_q)
  exact_lpd <- mll_outer + log(mean(exp(log_q - mll_outer)))

  list(elpd_i  = as.numeric(exact_lpd),
       se      = as.numeric(sd(log_q) / sqrt(D)),
       n_outer = D,
       n_inner = n_inner)
}

run_hsroc_pipeline <- function(rho, df_base,
                                rho_index = 1L,
                                seed_master = SEED_MASTER) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Fitting HSROC pipeline at rho = %.2f\n", rho))
  cat(sprintf("========================================\n"))

  df_eff <- apply_deff(df_base, rho)
  N      <- nrow(df_eff)

  cat(sprintf("Effective sample sizes at rho = %.2f:\n", rho))
  cat(sprintf("%-30s %6s %6s %6s %8s %8s\n",
              "Study", "Tests", "Pts", "DEff", "n_pos_eff", "n_neg_eff"))
  for (i in seq_len(nrow(df_eff))) {
    cat(sprintf("%-30s %6d %6d %6.2f %8d %8d\n",
                df_eff$study[i], df_eff$n_tests[i], df_eff$n_patients[i],
                df_eff$deff[i], df_eff$n_pos_eff[i], df_eff$n_neg_eff[i]))
  }
  sens_check <- max(abs(df_eff$TP_eff / pmax(df_eff$n_pos_eff, 1) -
                         df_eff$sens))
  spec_check <- max(abs(df_eff$TN_eff / pmax(df_eff$n_neg_eff, 1) -
                         df_eff$spec))
  cat(sprintf("Max |sens deviation| from rounding: %.4f\n", sens_check))
  cat(sprintf("Max |spec deviation| from rounding: %.4f\n", spec_check))

  jags_data <- list(
    N           = N,
    TP          = df_eff$TP_eff,
    FP          = df_eff$FP_eff,
    n_pos       = df_eff$n_pos_eff,
    n_neg       = df_eff$n_neg_eff,
    threshold_c = df_eff$threshold_c,
    G           = G,
    t_grid_c    = t_grid_c
  )

  set.seed(seed_master + rho_index * 1000L)
  inits_main <- make_inits(mcmc_main$chains,
                            base_seed = seed_master + rho_index * 1000L)

  cat(sprintf("[rho = %.2f] Main model: adapting (%d) ...\n",
              rho, mcmc_main$adapt))
  jags_model <- jags.model(
    textConnection(hsroc_model_string),
    data     = jags_data,
    inits    = inits_main,
    n.chains = mcmc_main$chains,
    n.adapt  = mcmc_main$adapt,
    quiet    = TRUE
  )
  cat(sprintf("[rho = %.2f] Burn-in (%d) ...\n", rho, mcmc_main$burnin))
  update(jags_model, n.iter = mcmc_main$burnin, progress.bar = "none")

  cat(sprintf("[rho = %.2f] Sampling (%d, thin = %d) ...\n",
              rho, mcmc_main$iter, mcmc_main$thin))

  posterior_global <- coda.samples(
    jags_model,
    variable.names = c("mu_theta", "gamma_0", "gamma_1", "beta",
                       "sigma_theta", "sigma_u"),
    n.iter = mcmc_main$iter, thin = mcmc_main$thin,
    progress.bar = "none"
  )
  posterior_grid <- coda.samples(
    jags_model,
    variable.names = c("sens_grid", "fpr_grid", "spec_grid"),
    n.iter = mcmc_main$iter, thin = mcmc_main$thin,
    progress.bar = "none"
  )

  cat(sprintf("\n[rho = %.2f] Posterior summary:\n", rho))
  print(summary(posterior_global))
  cat(sprintf("\n[rho = %.2f] Gelman-Rubin diagnostics:\n", rho))
  print(gelman.diag(posterior_global, multivariate = FALSE))

  cat(sprintf("\n[rho = %.2f] Quadratic spline sensitivity check ...\n", rho))
  jags_data_quad <- jags_data
  jags_data_quad$threshold_c2 <- df_eff$threshold_c^2
  jags_data_quad$G            <- NULL
  jags_data_quad$t_grid_c     <- NULL

  inits_quad <- make_inits(mcmc_quad$chains,
                            base_seed = seed_master + rho_index * 1000L + 500L)

  quad_result <- tryCatch({
    qm <- jags.model(textConnection(hsroc_quad_string),
                     data = jags_data_quad,
                     inits = inits_quad,
                     n.chains = mcmc_quad$chains,
                     n.adapt = mcmc_quad$adapt,
                     quiet = TRUE)
    update(qm, n.iter = mcmc_quad$burnin, progress.bar = "none")
    qp <- coda.samples(qm,
      variable.names = c("gamma_0", "gamma_1", "gamma_2", "mu_theta", "beta",
                         "sigma_theta", "sigma_u"),
      n.iter = mcmc_quad$iter, thin = mcmc_quad$thin,
      progress.bar = "none")
    qd <- do.call(rbind, qp)

    g2_med <- median(qd[, "gamma_2"])
    g2_lo  <- as.numeric(quantile(qd[, "gamma_2"], 0.025))
    g2_hi  <- as.numeric(quantile(qd[, "gamma_2"], 0.975))

    dic_lin  <- tryCatch(dic.samples(jags_model, n.iter = 10000,
                                      progress.bar = "none"),
                          error = function(e) NULL)
    dic_q    <- tryCatch(dic.samples(qm, n.iter = 10000,
                                      progress.bar = "none"),
                          error = function(e) NULL)

    list(gamma_2_med = g2_med,
         gamma_2_lo  = g2_lo,
         gamma_2_hi  = g2_hi,
         dic_linear  = if (!is.null(dic_lin)) sum(dic_lin$deviance) +
                          sum(dic_lin$penalty) else NA_real_,
         dic_quad    = if (!is.null(dic_q))   sum(dic_q$deviance) +
                          sum(dic_q$penalty)   else NA_real_)
  }, error = function(e) {
    cat(sprintf("  Quadratic model failed: %s\n", e$message))
    list(gamma_2_med = NA, gamma_2_lo = NA, gamma_2_hi = NA,
         dic_linear  = NA, dic_quad   = NA)
  })

  cat(sprintf("  gamma_2: median = %.4f, 95%% CrI = [%.4f, %.4f]\n",
              quad_result$gamma_2_med, quad_result$gamma_2_lo,
              quad_result$gamma_2_hi))
  if (!is.na(quad_result$dic_linear) && !is.na(quad_result$dic_quad)) {
    cat(sprintf("  DIC: linear = %.1f | quadratic = %.1f | delta = %.1f\n",
                quad_result$dic_linear, quad_result$dic_quad,
                quad_result$dic_quad - quad_result$dic_linear))
  }

  post_global <- do.call(rbind, posterior_global)
  post_grid   <- do.call(rbind, posterior_grid)

  mu_theta_draws    <- post_global[, "mu_theta"]
  gamma_0_draws     <- post_global[, "gamma_0"]
  gamma_1_draws     <- post_global[, "gamma_1"]
  beta_draws        <- post_global[, "beta"]
  sigma_theta_draws <- post_global[, "sigma_theta"]
  sigma_u_draws     <- post_global[, "sigma_u"]
  n_draws           <- nrow(post_global)

  sens_cols <- paste0("sens_grid[", seq_len(G), "]")
  spec_cols <- paste0("spec_grid[", seq_len(G), "]")
  fpr_cols  <- paste0("fpr_grid[",  seq_len(G), "]")

  sens_mat <- post_grid[, sens_cols]
  spec_mat <- post_grid[, spec_cols]
  fpr_mat  <- post_grid[, fpr_cols]

  sens_mean <- colMeans(sens_mat); spec_mean <- colMeans(spec_mat)
  fpr_mean  <- colMeans(fpr_mat)
  sens_lo   <- apply(sens_mat, 2, quantile, 0.025)
  sens_hi   <- apply(sens_mat, 2, quantile, 0.975)
  spec_lo   <- apply(spec_mat, 2, quantile, 0.025)
  spec_hi   <- apply(spec_mat, 2, quantile, 0.975)

  cat(sprintf("\n[rho = %.2f] Sampling log_lik for PSIS-LOO ...\n", rho))
  posterior_loglik <- coda.samples(
    jags_model,
    variable.names = "log_lik",
    n.iter = mcmc_main$iter, thin = mcmc_main$thin,
    progress.bar = "none"
  )
  n_chains_main <- mcmc_main$chains
  iter_per_chain <- nrow(posterior_loglik[[1]])
  log_lik_array <- array(NA_real_, dim = c(iter_per_chain, n_chains_main, N))
  ll_cols <- paste0("log_lik[", seq_len(N), "]")
  for (ch in seq_len(n_chains_main)) {
    chain_mat <- as.matrix(posterior_loglik[[ch]])
    log_lik_array[, ch, ] <- chain_mat[, ll_cols]
  }
  log_lik_mat <- matrix(log_lik_array,
                        nrow = iter_per_chain * n_chains_main, ncol = N)

  r_eff <- loo::relative_eff(exp(log_lik_array))
  loo_result <- loo::loo(log_lik_mat, r_eff = r_eff)

  cat(sprintf("\n[rho = %.2f] PSIS-LOO summary:\n", rho))
  print(loo_result)

  pareto_k <- loo_result$diagnostics$pareto_k
  elpd_loo_pointwise <- as.numeric(loo_result$pointwise[, "elpd_loo"])
  p_loo_pointwise    <- as.numeric(loo_result$pointwise[, "p_loo"])
  loo_study          <- df_eff$study

  cat(sprintf("\n[rho = %.2f] Per-study Pareto-k diagnostics:\n", rho))
  cat(sprintf("%-30s %10s %10s %10s %s\n",
              "Study", "elpd_loo", "p_loo", "pareto_k", "flag"))
  for (i in seq_len(N)) {
    flag <- if (pareto_k[i] > 1.0)  "*** k>1.0 ***"
            else if (pareto_k[i] > 0.7) "** k>0.7 **"
            else if (pareto_k[i] > 0.5) "  k>0.5"
            else                            "  ok"
    cat(sprintf("%-30s %10.2f %10.2f %10.2f  %s\n",
                loo_study[i], elpd_loo_pointwise[i],
                p_loo_pointwise[i], pareto_k[i], flag))
  }

  high_k_idx <- which(pareto_k > 0.7)
  exact_lpd_per_study <- rep(NA_real_, N)
  exact_lpd_se        <- rep(NA_real_, N)

  if (length(high_k_idx) > 0L) {
    cat(sprintf("\n[rho = %.2f] Pareto-k > 0.7 for %d studies; ",
                rho, length(high_k_idx)))
    cat("running exact LOO refits with n_inner = 1000 ...\n")
    for (idx_pos in seq_along(high_k_idx)) {
      i <- high_k_idx[idx_pos]
      cat(sprintf("  Exact LOO %d/%d: '%s' (k = %.2f) ... ",
                  idx_pos, length(high_k_idx),
                  df_eff$study[i], pareto_k[i]))
      res <- tryCatch(
        exact_loo_marginal(i, df_eff, hsroc_model_string,
                            t_grid_c, G, mcmc_loo,
                            seed_master, rho_index),
        error = function(e) {
          cat(sprintf("FAILED (%s)\n", e$message)); NULL
        })
      if (!is.null(res)) {
        exact_lpd_per_study[i] <- res$elpd_i
        exact_lpd_se[i]        <- res$se
        cat(sprintf("elpd = %.2f (se = %.2f)\n", res$elpd_i, res$se))
      }
    }
  } else {
    cat(sprintf("\n[rho = %.2f] All Pareto-k <= 0.7; PSIS-LOO is reliable.\n",
                rho))
  }

  elpd_loo_pointwise_final <- elpd_loo_pointwise
  for (i in seq_len(N)) {
    if (!is.na(exact_lpd_per_study[i])) {
      elpd_loo_pointwise_final[i] <- exact_lpd_per_study[i]
    }
  }

  elpd_loo_total    <- sum(elpd_loo_pointwise_final)
  elpd_loo_total_se <- sqrt(N * var(elpd_loo_pointwise_final))
  p_loo_total       <- as.numeric(loo_result$estimates["p_loo", "Estimate"])
  p_loo_total_se    <- as.numeric(loo_result$estimates["p_loo", "SE"])
  looic             <- as.numeric(loo_result$estimates["looic", "Estimate"])
  looic_se          <- as.numeric(loo_result$estimates["looic", "SE"])
  pareto_k_max      <- max(pareto_k, na.rm = TRUE)
  n_high_k          <- length(high_k_idx)

  cat(sprintf("\n[rho = %.2f] Final elpd_loo (with exact substitutions for k > 0.7): %.2f (se = %.2f)\n",
              rho, elpd_loo_total, elpd_loo_total_se))
  cat(sprintf("[rho = %.2f] p_loo = %.2f (se = %.2f) | p_loo / N = %.3f\n",
              rho, p_loo_total, p_loo_total_se, p_loo_total / N))
  cat(sprintf("[rho = %.2f] max Pareto-k = %.2f | studies with k > 0.7 = %d\n",
              rho, pareto_k_max, n_high_k))

  insample_lpd <- numeric(N)
  for (i in seq_len(N)) {
    ll_i <- log_lik_mat[, i]
    mll  <- max(ll_i)
    insample_lpd[i] <- mll + log(mean(exp(ll_i - mll)))
  }
  insample_total <- sum(insample_lpd)

  lambda_wide <- seq(-10, 10, length.out = 500)
  auc_draws   <- numeric(n_draws)
  for (d in seq_len(n_draws)) {
    al_d <- mu_theta_draws[d] * exp(-beta_draws[d] / 2)
    ps_d <- lambda_wide * exp(beta_draws[d] / 2)
    sens_d <- plogis((al_d + ps_d) / 2)
    fpr_d  <- plogis((-al_d + ps_d) / 2)
    o <- order(fpr_d)
    fs <- fpr_d[o]; ss <- sens_d[o]
    auc_draws[d] <- sum(diff(fs) * (ss[-1] + ss[-length(ss)]) / 2)
  }
  auc_sroc <- median(auc_draws)
  auc_lo   <- as.numeric(quantile(auc_draws, 0.025))
  auc_hi   <- as.numeric(quantile(auc_draws, 0.975))

  youden_vals   <- sens_mean + spec_mean - 1
  best_idx      <- which.max(youden_vals)
  youden_thresh <- t_grid_raw[best_idx]
  youden_sens   <- sens_mean[best_idx]
  youden_spec   <- spec_mean[best_idx]
  youden_fpr    <- fpr_mean[best_idx]
  youden_J      <- youden_vals[best_idx]
  youden_J_drs  <- sens_mat[, best_idx] + spec_mat[, best_idx] - 1
  youden_J_lo   <- as.numeric(quantile(youden_J_drs, 0.025))
  youden_J_hi   <- as.numeric(quantile(youden_J_drs, 0.975))

  max_sens_per_draw <- apply(sens_mat, 1, max)
  sens_target_candidates <- c(0.95, 0.90, 0.85, 0.80)
  sens_target <- NA
  for (tgt in sens_target_candidates) {
    if (mean(max_sens_per_draw >= tgt) >= 0.50) { sens_target <- tgt; break }
  }
  if (is.na(sens_target)) {
    sens_target <- as.numeric(round(quantile(max_sens_per_draw, 0.05), 2))
  }
  anchor_thresh_draws <- rep(NA_real_, n_draws)
  anchor_spec_draws   <- rep(NA_real_, n_draws)
  anchor_fpr_draws    <- rep(NA_real_, n_draws)
  for (d in seq_len(n_draws)) {
    s_curve <- sens_mat[d, ]
    idx_tgt <- suppressWarnings(max(which(s_curve >= sens_target)))
    if (is.finite(idx_tgt)) {
      anchor_thresh_draws[d] <- t_grid_raw[idx_tgt]
      anchor_spec_draws[d]   <- spec_mat[d, idx_tgt]
      anchor_fpr_draws[d]    <- fpr_mat[d,  idx_tgt]
    }
  }
  valid <- !is.na(anchor_thresh_draws)
  anchor_thresh_med <- median(anchor_thresh_draws, na.rm = TRUE)
  anchor_thresh_lo  <- as.numeric(quantile(anchor_thresh_draws, 0.025, na.rm = TRUE))
  anchor_thresh_hi  <- as.numeric(quantile(anchor_thresh_draws, 0.975, na.rm = TRUE))
  anchor_spec_med   <- median(anchor_spec_draws,  na.rm = TRUE)
  anchor_spec_lo    <- as.numeric(quantile(anchor_spec_draws, 0.025, na.rm = TRUE))
  anchor_spec_hi    <- as.numeric(quantile(anchor_spec_draws, 0.975, na.rm = TRUE))
  anchor_fpr_med    <- median(anchor_fpr_draws,  na.rm = TRUE)
  anchor_sens_label <- sprintf("%.0f%% Sens", sens_target * 100)

  max_spec_per_draw <- apply(spec_mat, 1, max)
  spec_target_candidates <- c(0.95, 0.90, 0.85, 0.80, 0.75)
  spec_target <- NA
  for (tgt in spec_target_candidates) {
    if (mean(max_spec_per_draw >= tgt) >= 0.50) { spec_target <- tgt; break }
  }
  if (is.na(spec_target)) {
    spec_target <- as.numeric(round(quantile(max_spec_per_draw, 0.05), 2))
  }
  spec_anchor_thresh_draws <- rep(NA_real_, n_draws)
  spec_anchor_sens_draws   <- rep(NA_real_, n_draws)
  spec_anchor_fpr_draws    <- rep(NA_real_, n_draws)
  for (d in seq_len(n_draws)) {
    sp_curve <- spec_mat[d, ]
    idx_sp   <- suppressWarnings(min(which(sp_curve >= spec_target)))
    if (is.finite(idx_sp)) {
      spec_anchor_thresh_draws[d] <- t_grid_raw[idx_sp]
      spec_anchor_sens_draws[d]   <- sens_mat[d, idx_sp]
      spec_anchor_fpr_draws[d]    <- fpr_mat[d,  idx_sp]
    }
  }
  valid_sp <- !is.na(spec_anchor_thresh_draws)
  spec90_thresh_med <- median(spec_anchor_thresh_draws, na.rm = TRUE)
  spec90_thresh_lo  <- as.numeric(quantile(spec_anchor_thresh_draws, 0.025, na.rm = TRUE))
  spec90_thresh_hi  <- as.numeric(quantile(spec_anchor_thresh_draws, 0.975, na.rm = TRUE))
  spec90_sens_med   <- median(spec_anchor_sens_draws,  na.rm = TRUE)
  spec90_sens_lo    <- as.numeric(quantile(spec_anchor_sens_draws, 0.025, na.rm = TRUE))
  spec90_sens_hi    <- as.numeric(quantile(spec_anchor_sens_draws, 0.975, na.rm = TRUE))
  spec90_fpr_med    <- median(spec_anchor_fpr_draws,   na.rm = TRUE)
  anchor_spec_label <- sprintf("%.0f%% Spec", spec_target * 100)

  logit_sens_summary_draws <- (mu_theta_draws * exp(-beta_draws / 2) +
                                gamma_0_draws  * exp( beta_draws / 2)) / 2
  logit_fpr_summary_draws  <- (-mu_theta_draws * exp(-beta_draws / 2) +
                                 gamma_0_draws  * exp( beta_draws / 2)) / 2
  sens_summary_draws <- plogis(logit_sens_summary_draws)
  fpr_summary_draws  <- plogis(logit_fpr_summary_draws)

  summary_sens   <- median(sens_summary_draws)
  summary_fpr    <- median(fpr_summary_draws)
  summary_spec   <- 1 - summary_fpr
  summary_sens_lo <- as.numeric(quantile(sens_summary_draws, 0.025))
  summary_sens_hi <- as.numeric(quantile(sens_summary_draws, 0.975))
  summary_fpr_lo  <- as.numeric(quantile(fpr_summary_draws, 0.025))
  summary_fpr_hi  <- as.numeric(quantile(fpr_summary_draws, 0.975))
  summary_spec_lo <- as.numeric(quantile(1 - fpr_summary_draws, 0.025))
  summary_spec_hi <- as.numeric(quantile(1 - fpr_summary_draws, 0.975))

  sens_youden_vec <- numeric(n_draws); spec_youden_vec <- numeric(n_draws)
  sens_low_vec    <- rep(NA_real_, n_draws)
  spec_low_vec    <- rep(NA_real_, n_draws)
  sens_high_vec   <- rep(NA_real_, n_draws)
  spec_high_vec   <- rep(NA_real_, n_draws)
  for (d in seq_len(n_draws)) {
    j_d  <- sens_mat[d, ] + spec_mat[d, ] - 1
    bd   <- which.max(j_d)
    sens_youden_vec[d] <- sens_mat[d, bd]
    spec_youden_vec[d] <- spec_mat[d, bd]

    idx_lo <- suppressWarnings(max(which(sens_mat[d, ] >= sens_target)))
    if (is.finite(idx_lo)) {
      sens_low_vec[d] <- sens_mat[d, idx_lo]
      spec_low_vec[d] <- spec_mat[d, idx_lo]
    }
    idx_hi <- suppressWarnings(min(which(spec_mat[d, ] >= spec_target)))
    if (is.finite(idx_hi)) {
      sens_high_vec[d] <- sens_mat[d, idx_hi]
      spec_high_vec[d] <- spec_mat[d, idx_hi]
    }
  }

  out <- list(
    rho = rho, N = N, df_eff = df_eff,

    mu_theta_draws = mu_theta_draws, gamma_0_draws = gamma_0_draws,
    gamma_1_draws  = gamma_1_draws,  beta_draws    = beta_draws,
    sigma_theta_draws = sigma_theta_draws, sigma_u_draws = sigma_u_draws,
    n_draws = n_draws,

    sens_mat = sens_mat, spec_mat = spec_mat, fpr_mat = fpr_mat,
    sens_mean = sens_mean, spec_mean = spec_mean, fpr_mean = fpr_mean,
    sens_lo = sens_lo, sens_hi = sens_hi,
    spec_lo = spec_lo, spec_hi = spec_hi,

    auc_sroc = auc_sroc, auc_lo = auc_lo, auc_hi = auc_hi,
    auc_draws = auc_draws,

    youden_thresh = youden_thresh, youden_sens = youden_sens,
    youden_spec   = youden_spec,   youden_fpr  = youden_fpr,
    youden_J = youden_J, youden_J_lo = youden_J_lo, youden_J_hi = youden_J_hi,
    best_idx = best_idx,

    sens_target = sens_target, anchor_sens_label = anchor_sens_label,
    anchor_thresh_med = anchor_thresh_med,
    anchor_thresh_lo  = anchor_thresh_lo,
    anchor_thresh_hi  = anchor_thresh_hi,
    anchor_spec_med = anchor_spec_med,
    anchor_spec_lo  = anchor_spec_lo,
    anchor_spec_hi  = anchor_spec_hi,
    anchor_fpr_med  = anchor_fpr_med,
    anchor_thresh_draws = anchor_thresh_draws,
    anchor_spec_draws   = anchor_spec_draws,
    anchor_fpr_draws    = anchor_fpr_draws,
    valid_anchor_sens   = valid,

    spec_target = spec_target, anchor_spec_label = anchor_spec_label,
    spec90_thresh_med = spec90_thresh_med,
    spec90_thresh_lo  = spec90_thresh_lo,
    spec90_thresh_hi  = spec90_thresh_hi,
    spec90_sens_med = spec90_sens_med,
    spec90_sens_lo  = spec90_sens_lo,
    spec90_sens_hi  = spec90_sens_hi,
    spec90_fpr_med  = spec90_fpr_med,
    spec_anchor_thresh_draws = spec_anchor_thresh_draws,
    spec_anchor_sens_draws   = spec_anchor_sens_draws,
    spec_anchor_fpr_draws    = spec_anchor_fpr_draws,
    valid_anchor_spec        = valid_sp,

    summary_sens = summary_sens, summary_spec = summary_spec,
    summary_fpr  = summary_fpr,
    summary_sens_lo = summary_sens_lo, summary_sens_hi = summary_sens_hi,
    summary_fpr_lo  = summary_fpr_lo,  summary_fpr_hi  = summary_fpr_hi,
    summary_spec_lo = summary_spec_lo, summary_spec_hi = summary_spec_hi,
    sens_summary_draws = sens_summary_draws,
    fpr_summary_draws  = fpr_summary_draws,
    logit_sens_summary_draws = logit_sens_summary_draws,
    logit_fpr_summary_draws  = logit_fpr_summary_draws,

    sens_youden_vec = sens_youden_vec, spec_youden_vec = spec_youden_vec,
    sens_low_vec    = sens_low_vec,    spec_low_vec    = spec_low_vec,
    sens_high_vec   = sens_high_vec,   spec_high_vec   = spec_high_vec,

    quad_result = quad_result,

    loo_study                 = loo_study,
    loo_result                = loo_result,
    pareto_k                  = pareto_k,
    elpd_loo_pointwise        = elpd_loo_pointwise,
    p_loo_pointwise           = p_loo_pointwise,
    high_k_idx                = high_k_idx,
    exact_lpd_per_study       = exact_lpd_per_study,
    exact_lpd_se              = exact_lpd_se,
    elpd_loo_pointwise_final  = elpd_loo_pointwise_final,
    elpd_loo_total            = elpd_loo_total,
    elpd_loo_total_se         = elpd_loo_total_se,
    p_loo_total               = p_loo_total,
    p_loo_total_se            = p_loo_total_se,
    looic                     = looic,
    looic_se                  = looic_se,
    pareto_k_max              = pareto_k_max,
    n_high_k                  = n_high_k,
    insample_lpd              = insample_lpd,
    insample_total            = insample_total
  )

  return(out)
}

results <- list()
for (k in seq_along(rho_values)) {
  rho_k <- rho_values[k]
  results[[sprintf("%.2f", rho_k)]] <- run_hsroc_pipeline(
    rho          = rho_k,
    df_base      = df_raw,
    rho_index    = k,
    seed_master  = SEED_MASTER
  )
}

for (k in seq_along(rho_values)) {
  r <- results[[k]]
  study_grid_idx <- vapply(r$df_eff$threshold,
                            function(t) which.min(abs(t_grid_raw - t)),
                            integer(1))
  sens_marg_draws <- rowMeans(r$sens_mat[, study_grid_idx])
  spec_marg_draws <- rowMeans(r$spec_mat[, study_grid_idx])
  fpr_marg_draws  <- 1 - spec_marg_draws

  r$marg_summary_sens     <- median(sens_marg_draws)
  r$marg_summary_spec     <- median(spec_marg_draws)
  r$marg_summary_fpr      <- median(fpr_marg_draws)
  r$marg_summary_sens_lo  <- as.numeric(quantile(sens_marg_draws, 0.025))
  r$marg_summary_sens_hi  <- as.numeric(quantile(sens_marg_draws, 0.975))
  r$marg_summary_spec_lo  <- as.numeric(quantile(spec_marg_draws, 0.025))
  r$marg_summary_spec_hi  <- as.numeric(quantile(spec_marg_draws, 0.975))
  r$marg_summary_fpr_lo   <- as.numeric(quantile(fpr_marg_draws, 0.025))
  r$marg_summary_fpr_hi   <- as.numeric(quantile(fpr_marg_draws, 0.975))
  r$sens_marg_draws       <- sens_marg_draws
  r$spec_marg_draws       <- spec_marg_draws
  r$fpr_marg_draws        <- fpr_marg_draws

  r$p_loo_blended <- r$insample_total - r$elpd_loo_total

  results[[k]] <- r
}

fmt_ci <- function(med, lo, hi, dec = 3) {
  sprintf(paste0("%.", dec, "f [%.", dec, "f-%.", dec, "f]"), med, lo, hi)
}
fmt_pct <- function(x, dec = 1) sprintf(paste0("%.", dec, "f%%"), x * 100)

cat("\n============================================================\n")
cat("rho SENSITIVITY SWEEP - SUMMARY TABLE\n")
cat("============================================================\n")
cat(sprintf("Primary analysis: rho = %.2f\n\n", primary_rho))

comp <- data.frame(
  rho                  = rho_values,
  primary              = ifelse(rho_values == primary_rho, "*", " "),
  mean_DEff            = sapply(results, function(r) round(mean(r$df_eff$deff), 2)),
  total_n_eff          = sapply(results, function(r) sum(r$df_eff$n_pos_eff +
                                                          r$df_eff$n_neg_eff)),
  AUC                  = sapply(results, function(r) fmt_ci(r$auc_sroc,
                                                             r$auc_lo, r$auc_hi)),
  Youden_cutoff_pct    = sapply(results, function(r) sprintf("%.3f%%",
                                                              r$youden_thresh * 100)),
  Youden_sens          = sapply(results, function(r) sprintf("%.3f", r$youden_sens)),
  Youden_spec          = sapply(results, function(r) sprintf("%.3f", r$youden_spec)),
  Youden_J             = sapply(results, function(r) fmt_ci(r$youden_J,
                                                             r$youden_J_lo,
                                                             r$youden_J_hi)),
  sens_target          = sapply(results, function(r) fmt_pct(r$sens_target, 0)),
  sens_anchor_cutoff   = sapply(results, function(r) sprintf("%.3f%% [%.3f-%.3f]",
                                                              r$anchor_thresh_med * 100,
                                                              r$anchor_thresh_lo * 100,
                                                              r$anchor_thresh_hi * 100)),
  spec_at_sens_anchor  = sapply(results, function(r) fmt_ci(r$anchor_spec_med,
                                                             r$anchor_spec_lo,
                                                             r$anchor_spec_hi)),
  spec_target          = sapply(results, function(r) fmt_pct(r$spec_target, 0)),
  spec_anchor_cutoff   = sapply(results, function(r) sprintf("%.3f%% [%.3f-%.3f]",
                                                              r$spec90_thresh_med * 100,
                                                              r$spec90_thresh_lo * 100,
                                                              r$spec90_thresh_hi * 100)),
  sens_at_spec_anchor  = sapply(results, function(r) fmt_ci(r$spec90_sens_med,
                                                             r$spec90_sens_lo,
                                                             r$spec90_sens_hi)),

  op_pt_gamma0_sens    = sapply(results, function(r) fmt_ci(r$summary_sens,
                                                             r$summary_sens_lo,
                                                             r$summary_sens_hi)),
  op_pt_gamma0_spec    = sapply(results, function(r) fmt_ci(r$summary_spec,
                                                             r$summary_spec_lo,
                                                             r$summary_spec_hi)),

  summary_marg_sens    = sapply(results, function(r) fmt_ci(r$marg_summary_sens,
                                                             r$marg_summary_sens_lo,
                                                             r$marg_summary_sens_hi)),
  summary_marg_spec    = sapply(results, function(r) fmt_ci(r$marg_summary_spec,
                                                             r$marg_summary_spec_lo,
                                                             r$marg_summary_spec_hi)),
  elpd_loo             = sapply(results, function(r) sprintf("%.2f (%.2f)",
                                                                r$elpd_loo_total,
                                                                r$elpd_loo_total_se)),
  p_loo_psis           = sapply(results, function(r) sprintf("%.2f (%.2f)",
                                                                r$p_loo_total,
                                                                r$p_loo_total_se)),
  p_loo_blended        = sapply(results, function(r) sprintf("%.2f",
                                                                r$p_loo_blended)),
  looic                = sapply(results, function(r) sprintf("%.2f (%.2f)",
                                                                r$looic, r$looic_se)),
  pareto_k_max         = sapply(results, function(r) sprintf("%.2f", r$pareto_k_max)),
  n_high_k             = sapply(results, function(r) r$n_high_k),
  insample_lpd_total   = sapply(results, function(r) sprintf("%.2f", r$insample_total)),
  stringsAsFactors     = FALSE
)
print(comp, row.names = FALSE)

primary <- results[[sprintf("%.2f", primary_rho)]]

df              <- primary$df_eff
N               <- primary$N

mu_theta_draws  <- primary$mu_theta_draws
gamma_0_draws   <- primary$gamma_0_draws
gamma_1_draws   <- primary$gamma_1_draws
beta_draws      <- primary$beta_draws
sigma_theta_draws <- primary$sigma_theta_draws
sigma_u_draws     <- primary$sigma_u_draws
n_draws         <- primary$n_draws
post_global     <- cbind(mu_theta = mu_theta_draws,
                          gamma_0  = gamma_0_draws,
                          gamma_1  = gamma_1_draws,
                          beta     = beta_draws,
                          sigma_theta = sigma_theta_draws,
                          sigma_u     = sigma_u_draws)

sens_mat <- primary$sens_mat; spec_mat <- primary$spec_mat
fpr_mat  <- primary$fpr_mat
sens_mean <- primary$sens_mean; spec_mean <- primary$spec_mean
fpr_mean  <- primary$fpr_mean
sens_lo <- primary$sens_lo; sens_hi <- primary$sens_hi
spec_lo <- primary$spec_lo; spec_hi <- primary$spec_hi

auc_sroc <- primary$auc_sroc; auc_lo <- primary$auc_lo; auc_hi <- primary$auc_hi
auc_draws <- primary$auc_draws

youden_thresh <- primary$youden_thresh
youden_sens   <- primary$youden_sens
youden_spec   <- primary$youden_spec
youden_fpr    <- primary$youden_fpr
youden_J      <- primary$youden_J
youden_J_lo   <- primary$youden_J_lo
youden_J_hi   <- primary$youden_J_hi
best_idx      <- primary$best_idx

sens_target          <- primary$sens_target
anchor_sens_target   <- sens_target
anchor_sens_label    <- primary$anchor_sens_label
anchor_thresh_med    <- primary$anchor_thresh_med
anchor_thresh_lo     <- primary$anchor_thresh_lo
anchor_thresh_hi     <- primary$anchor_thresh_hi
anchor_spec_med      <- primary$anchor_spec_med
anchor_spec_lo       <- primary$anchor_spec_lo
anchor_spec_hi       <- primary$anchor_spec_hi
anchor_fpr_med       <- primary$anchor_fpr_med
anchor_thresh_draws  <- primary$anchor_thresh_draws
anchor_spec_draws    <- primary$anchor_spec_draws
anchor_fpr_draws     <- primary$anchor_fpr_draws
valid                <- primary$valid_anchor_sens

spec_target              <- primary$spec_target
anchor_spec_target       <- spec_target
anchor_spec_label        <- primary$anchor_spec_label
spec90_thresh_med        <- primary$spec90_thresh_med
spec90_thresh_lo         <- primary$spec90_thresh_lo
spec90_thresh_hi         <- primary$spec90_thresh_hi
spec90_sens_med          <- primary$spec90_sens_med
spec90_sens_lo           <- primary$spec90_sens_lo
spec90_sens_hi           <- primary$spec90_sens_hi
spec90_fpr_med           <- primary$spec90_fpr_med
spec_anchor_thresh_draws <- primary$spec_anchor_thresh_draws
spec_anchor_sens_draws   <- primary$spec_anchor_sens_draws
spec_anchor_fpr_draws    <- primary$spec_anchor_fpr_draws
valid_sp                 <- primary$valid_anchor_spec

gamma0_sens     <- primary$summary_sens
gamma0_spec     <- primary$summary_spec
gamma0_fpr      <- primary$summary_fpr
gamma0_sens_lo  <- primary$summary_sens_lo
gamma0_sens_hi  <- primary$summary_sens_hi
gamma0_fpr_lo   <- primary$summary_fpr_lo
gamma0_fpr_hi   <- primary$summary_fpr_hi
gamma0_spec_lo  <- primary$summary_spec_lo
gamma0_spec_hi  <- primary$summary_spec_hi
sens_gamma0_draws       <- primary$sens_summary_draws
fpr_gamma0_draws        <- primary$fpr_summary_draws
logit_sens_gamma0_draws <- primary$logit_sens_summary_draws
logit_fpr_gamma0_draws  <- primary$logit_fpr_summary_draws

summary_sens     <- primary$marg_summary_sens
summary_spec     <- primary$marg_summary_spec
summary_fpr      <- primary$marg_summary_fpr
summary_sens_lo  <- primary$marg_summary_sens_lo
summary_sens_hi  <- primary$marg_summary_sens_hi
summary_spec_lo  <- primary$marg_summary_spec_lo
summary_spec_hi  <- primary$marg_summary_spec_hi
summary_fpr_lo   <- primary$marg_summary_fpr_lo
summary_fpr_hi   <- primary$marg_summary_fpr_hi
sens_summary_draws <- primary$sens_marg_draws
fpr_summary_draws  <- primary$fpr_marg_draws
spec_summary_draws <- primary$spec_marg_draws

sens_youden_vec <- primary$sens_youden_vec
spec_youden_vec <- primary$spec_youden_vec
sens_low_vec    <- primary$sens_low_vec
spec_low_vec    <- primary$spec_low_vec
sens_high_vec   <- primary$sens_high_vec
spec_high_vec   <- primary$spec_high_vec

loo_study                <- primary$loo_study
loo_result               <- primary$loo_result
pareto_k                 <- primary$pareto_k
elpd_loo_pointwise       <- primary$elpd_loo_pointwise
p_loo_pointwise          <- primary$p_loo_pointwise
elpd_loo_pointwise_final <- primary$elpd_loo_pointwise_final
exact_lpd_per_study      <- primary$exact_lpd_per_study
elpd_loo_total           <- primary$elpd_loo_total
elpd_loo_total_se        <- primary$elpd_loo_total_se
p_loo_total              <- primary$p_loo_total
p_loo_total_se           <- primary$p_loo_total_se
looic                    <- primary$looic
looic_se                 <- primary$looic_se
pareto_k_max             <- primary$pareto_k_max
n_high_k                 <- primary$n_high_k
insample_lpd             <- primary$insample_lpd
insample_total           <- primary$insample_total
p_loo_blended            <- primary$p_loo_blended

cat(sprintf("\nPrimary results loaded for figures: rho = %.2f\n", primary_rho))
cat(sprintf("  AUC = %.3f [%.3f-%.3f]\n", auc_sroc, auc_lo, auc_hi))
cat(sprintf("  Youden cutoff = %.3f%% (sens = %.3f, spec = %.3f)\n",
            youden_thresh * 100, youden_sens, youden_spec))
cat(sprintf("  Sens target = %.0f%%; cutoff = %.3f%% (CrI %.3f-%.3f); spec at anchor = %.3f\n",
            sens_target * 100, anchor_thresh_med * 100,
            anchor_thresh_lo * 100, anchor_thresh_hi * 100, anchor_spec_med))
cat(sprintf("  Spec target = %.0f%%; cutoff = %.3f%% (CrI %.3f-%.3f); sens at anchor = %.3f\n",
            spec_target * 100, spec90_thresh_med * 100,
            spec90_thresh_lo * 100, spec90_thresh_hi * 100, spec90_sens_med))

finite_marg <- is.finite(sens_summary_draws) & is.finite(fpr_summary_draws)
conf_draws_df <- data.frame(
  fpr  = fpr_summary_draws[finite_marg],
  sens = sens_summary_draws[finite_marg]
)

set.seed(SEED_MASTER + 7L)
n_pred     <- length(mu_theta_draws)
theta_new  <- rnorm(n_pred, mu_theta_draws, sigma_theta_draws)
u_new      <- rnorm(n_pred, 0, sigma_u_draws)
lambda_new <- gamma_0_draws + u_new
logit_sens_pred <- (theta_new  * exp(-beta_draws / 2) +
                     lambda_new * exp( beta_draws / 2)) / 2
logit_fpr_pred  <- (-theta_new  * exp(-beta_draws / 2) +
                      lambda_new * exp( beta_draws / 2)) / 2
finite_pred <- is.finite(logit_fpr_pred) & is.finite(logit_sens_pred)
pred_draws_df <- data.frame(
  fpr  = plogis(logit_fpr_pred[finite_pred]),
  sens = plogis(logit_sens_pred[finite_pred])
)

set.seed(SEED_MASTER + 8L)
conf_thin <- if (nrow(conf_draws_df) > 5000)
  conf_draws_df[sample.int(nrow(conf_draws_df), 5000), ] else conf_draws_df
pred_thin <- if (nrow(pred_draws_df) > 5000)
  pred_draws_df[sample.int(nrow(pred_draws_df), 5000), ] else pred_draws_df

sroc_df <- data.frame(
  fpr = fpr_mean, sens = sens_mean,
  sens_lo = sens_lo, sens_hi = sens_hi
) %>% arrange(fpr)

study_pts <- data.frame(
  fpr = df$fpr, sens = df$sens, study = df$study,
  n   = df$n_pos + df$n_neg
)

param_label <- sprintf(
  paste0("HSROC parameters (posterior medians) | rho_cluster = %.2f\n",
         "mu_theta = %.3f  gamma_0 = %.3f  gamma_1 = %.3f\n",
         "beta = %.3f  sigma_theta = %.3f  sigma_u = %.3f\n",
         "AUC = %.3f [%.3f-%.3f]"),
  primary_rho,
  median(mu_theta_draws), median(gamma_0_draws), median(gamma_1_draws),
  median(beta_draws), median(sigma_theta_draws), median(sigma_u_draws),
  auc_sroc, auc_lo, auc_hi)

hsroc_label_list <- list(
  data.frame(x = summary_fpr, y = summary_sens,
             label = sprintf("Summary point (marginalised)\nSens = %.2f, Spec = %.2f",
                             summary_sens, summary_spec),
             color = "#D94F3D", stringsAsFactors = FALSE),
  data.frame(x = gamma0_fpr, y = gamma0_sens,
             label = sprintf("Operating point at mean observed\nthreshold (γ₀)\nSens = %.2f, Spec = %.2f",
                             gamma0_sens, gamma0_spec),
             color = "#9C2A1F", stringsAsFactors = FALSE),
  data.frame(x = youden_fpr, y = youden_sens,
             label = sprintf("Youden (J = %.2f)\nCutoff = %.2f%%\nSens = %.2f, Spec = %.2f",
                             youden_J, youden_thresh * 100, youden_sens, youden_spec),
             color = "#E8850C", stringsAsFactors = FALSE)
)
if (!is.na(anchor_fpr_med)) {
  hsroc_label_list[[length(hsroc_label_list) + 1]] <- data.frame(
    x = anchor_fpr_med, y = sens_target,
    label = sprintf("%s Anchor\nCutoff = %.2f%%\nSpec = %.2f [%.2f-%.2f]",
                    anchor_sens_label, anchor_thresh_med * 100,
                    anchor_spec_med, anchor_spec_lo, anchor_spec_hi),
    color = "#2A9D8F", stringsAsFactors = FALSE)
}
if (!is.na(spec90_fpr_med)) {
  hsroc_label_list[[length(hsroc_label_list) + 1]] <- data.frame(
    x = spec90_fpr_med, y = spec90_sens_med,
    label = sprintf("%s Anchor\nCutoff = %.2f%%\nSens = %.2f [%.2f-%.2f]",
                    anchor_spec_label, spec90_thresh_med * 100,
                    spec90_sens_med, spec90_sens_lo, spec90_sens_hi),
    color = "#7B2D8E", stringsAsFactors = FALSE)
}
hsroc_labels <- do.call(rbind, hsroc_label_list)

fig1 <- ggplot() +
  stat_density_2d(data = pred_thin, aes(x = fpr, y = sens),
                  geom = "polygon", bins = 4,
                  fill = "#1B4F8A", alpha = 0.04,
                  color = "#1B4F8A", linetype = "longdash", linewidth = 0.4) +
  stat_density_2d(data = conf_thin, aes(x = fpr, y = sens),
                  geom = "polygon", bins = 6,
                  fill = "#1B4F8A", alpha = 0.10,
                  color = "#1B4F8A", linetype = "solid", linewidth = 0.5) +
  geom_ribbon(data = sroc_df, aes(x = fpr, ymin = sens_lo, ymax = sens_hi),
              fill = "#1B4F8A", alpha = 0.10) +
  geom_line(data = sroc_df, aes(x = fpr, y = sens),
            color = "#1B4F8A", linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.5) +
  geom_point(data = study_pts, aes(x = fpr, y = sens, size = n),
             color = "#444444", alpha = 0.70, shape = 16) +
  scale_size_continuous(range = c(2, 6), name = "Sample size",
                        guide = guide_legend(override.aes = list(alpha = 0.8))) +

  geom_point(aes(x = summary_fpr, y = summary_sens),
             color = "#D94F3D", fill = "#D94F3D", size = 5, shape = 18) +

  geom_point(aes(x = gamma0_fpr, y = gamma0_sens),
             color = "#9C2A1F", fill = "#9C2A1F",
             size = 2.6, shape = 21, stroke = 0.6) +
  geom_point(aes(x = youden_fpr, y = youden_sens),
             color = "#E8850C", fill = "#E8850C", size = 4, shape = 17) +
  geom_point(aes(x = anchor_fpr_med, y = sens_target),
             color = "#2A9D8F", fill = "#2A9D8F", size = 4, shape = 23) +
  geom_point(aes(x = spec90_fpr_med, y = spec90_sens_med),
             color = "#7B2D8E", fill = "#7B2D8E", size = 4, shape = 15) +
  geom_label_repel(data = hsroc_labels,
                   aes(x = x, y = y, label = label),
                   color = hsroc_labels$color, size = 2.6,
                   fontface = "bold", lineheight = 0.9,
                   fill = alpha("white", 0.85), label.size = 0.2,
                   box.padding = 0.6, point.padding = 0.4,
                   segment.color = "grey50", segment.size = 0.3,
                   min.segment.length = 0.2, max.overlaps = 20,
                   seed = SEED_MASTER) +
  annotate("text", x = 0.72, y = 0.22,
           label = "----  95% Credible region\n- - -  95% Prediction region",
           color = "#1B4F8A", size = 2.8, hjust = 0, lineheight = 1.5) +
  geom_label(data = data.frame(x = 0.02, y = 0.10, lab = param_label),
             aes(x = x, y = y, label = lab),
             hjust = 0, vjust = 0, size = 2.5, fill = "#EEF4FB",
             color = "#1B4F8A", label.size = 0.4, lineheight = 1.3,
             inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  labs(
    title    = "HSROC - dd-cfDNA for Cardiac Allograft Rejection",
    subtitle = sprintf(
      "Bayesian threshold-regression HSROC | %d studies | rho_cluster = %.2f (primary) | AUC = %.3f [%.3f-%.3f]",
      N, primary_rho, auc_sroc, auc_lo, auc_hi),
    x = "False Positive Rate (1 - Specificity)",
    y = "Sensitivity (True Positive Rate)",
    caption = "rho_cluster sensitivity reported in rho_sensitivity_summary.csv"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.background  = element_rect(fill = "#FAFAFA"),
        plot.title        = element_text(face = "bold", size = 12),
        plot.subtitle     = element_text(size = 9, color = "grey40"),
        plot.caption      = element_text(size = 7.5, color = "grey50", hjust = 0),
        legend.position   = c(0.85, 0.15),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text       = element_text(size = 8),
        legend.title      = element_text(size = 9, face = "bold"))

print(fig1)
cat("Figure 1 (HSROC) saved.\n")

thresh_df <- data.frame(
  thresh = t_grid_raw,
  sens = sens_mean, sens_lo = sens_lo, sens_hi = sens_hi,
  spec = spec_mean, spec_lo = spec_lo, spec_hi = spec_hi
) %>%
  pivot_longer(cols = c(sens, spec),
               names_to = "metric", values_to = "value") %>%
  mutate(
    lo     = ifelse(metric == "sens", sens_lo, spec_lo),
    hi     = ifelse(metric == "sens", sens_hi, spec_hi),
    metric = ifelse(metric == "sens", "Sensitivity", "Specificity")
  )

nudge_1b <- (max(t_grid_raw) - min(t_grid_raw)) * 100 * 0.06

fig1b <- ggplot(thresh_df, aes(x = thresh * 100, y = value,
                                color = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = youden_thresh * 100,
             color = "#D94F3D", linetype = "dashed", linewidth = 0.9) +
  geom_label(data = data.frame(
               x   = youden_thresh * 100 + nudge_1b,
               y   = 0.12,
               lab = sprintf("Youden\n%.2f%%", youden_thresh * 100)),
             aes(x = x, y = y, label = lab),
             color = "#D94F3D", fill = alpha("white", 0.85),
             size = 3, hjust = 0, label.size = 0.2, inherit.aes = FALSE)

if (!is.na(anchor_thresh_med)) {
  fig1b <- fig1b +
    geom_vline(xintercept = anchor_thresh_med * 100,
               color = "#2A9D8F", linetype = "dashed", linewidth = 0.9) +
    geom_hline(yintercept = sens_target, linetype = "dotted",
               color = "#2A9D8F", linewidth = 0.6) +
    geom_label(data = data.frame(
                 x   = anchor_thresh_med * 100 + nudge_1b,
                 y   = 0.25,
                 lab = sprintf("%s\n%.2f%%", anchor_sens_label,
                               anchor_thresh_med * 100)),
               aes(x = x, y = y, label = lab),
               color = "#2A9D8F", fill = alpha("white", 0.85),
               size = 3, hjust = 0, label.size = 0.2,
               inherit.aes = FALSE)
}
if (!is.na(spec90_thresh_med)) {
  fig1b <- fig1b +
    geom_vline(xintercept = spec90_thresh_med * 100,
               color = "#7B2D8E", linetype = "dashed", linewidth = 0.9) +
    geom_hline(yintercept = spec_target, linetype = "dotted",
               color = "#7B2D8E", linewidth = 0.6) +
    geom_label(data = data.frame(
                 x   = spec90_thresh_med * 100 + nudge_1b,
                 y   = 0.38,
                 lab = sprintf("%s\n%.2f%%", anchor_spec_label,
                               spec90_thresh_med * 100)),
               aes(x = x, y = y, label = lab),
               color = "#7B2D8E", fill = alpha("white", 0.85),
               size = 3, hjust = 0, label.size = 0.2,
               inherit.aes = FALSE)
}

fig1b <- fig1b +
  scale_color_manual(values = c("Sensitivity" = "#1B4F8A",
                                "Specificity" = "#2A9D8F")) +
  scale_fill_manual(values  = c("Sensitivity" = "#1B4F8A",
                                "Specificity" = "#2A9D8F")) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title    = "Sensitivity & Specificity vs dd-cfDNA Threshold",
       subtitle = sprintf("Posterior means with 95%% CrI | rho_cluster = %.2f (primary)",
                          primary_rho),
       x = "dd-cfDNA Threshold (%)", y = "Probability",
       color = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.background = element_rect(fill = "#FAFAFA"),
        legend.position  = "top",
        plot.title       = element_text(face = "bold", size = 11),
        plot.subtitle    = element_text(size = 9, color = "grey40"))

print(fig1b)
cat("Figure 1b saved.\n")

jeffreys_ci <- function(k, n, level = 0.95) {
  alpha <- (1 - level) / 2
  list(lo = qbeta(alpha, k + 0.5, n - k + 0.5),
       hi = qbeta(1 - alpha, k + 0.5, n - k + 0.5))
}

forest_data <- df %>%
  mutate(
    sens_lo = mapply(function(k, n) jeffreys_ci(k, n)$lo, TP, n_pos),
    sens_hi = mapply(function(k, n) jeffreys_ci(k, n)$hi, TP, n_pos),
    spec_lo = mapply(function(k, n) jeffreys_ci(k, n)$lo, TN, n_neg),
    spec_hi = mapply(function(k, n) jeffreys_ci(k, n)$hi, TN, n_neg),
    pt_size_sens = 2 + 6 * ((1 / (sens_hi - sens_lo + 0.01)) /
                              max(1 / (sens_hi - sens_lo + 0.01))),
    pt_size_spec = 2 + 6 * ((1 / (spec_hi - spec_lo + 0.01)) /
                              max(1 / (spec_hi - spec_lo + 0.01))),
    label_sens   = sprintf("%.3f [%.3f-%.3f]", sens, sens_lo, sens_hi),
    label_spec   = sprintf("%.3f [%.3f-%.3f]", spec, spec_lo, spec_hi),
    study_f      = factor(study, levels = unique(rev(study)))
  )

make_forest <- function(data, metric, lo, hi, pt_size, val_label,
                        xlab, title, color, note) {
  ggplot(data, aes(y = study_f)) +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +

    geom_errorbar(aes(xmin = .data[[lo]], xmax = .data[[hi]]),
                  width = 0.25, color = color, linewidth = 0.8,
                  orientation = "y") +
    geom_point(aes(x = .data[[metric]], size = .data[[pt_size]]),
               color = color, alpha = 0.85) +
    geom_text(aes(x = 1.15, label = .data[[val_label]]),
              size = 2.9, hjust = 0, color = "grey20") +
    scale_size_identity() +
    scale_x_continuous(limits = c(-0.05, 1.55),
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    annotate("text", x = 1.15, y = nrow(data) + 0.8,
             label = "Value [95% CI]", hjust = 0, size = 3, fontface = "bold") +
    annotate("text", x = -0.05, y = -0.5, label = note,
             hjust = 0, size = 2.5, color = "grey50", fontface = "italic") +
    labs(title = title, x = xlab, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.background   = element_rect(fill = "#FAFAFA"),
          panel.grid.major.y = element_blank(),
          legend.position    = "none",
          plot.title         = element_text(face = "bold", size = 10),
          axis.text.y        = element_text(size = 8.5))
}

fig2 <- make_forest(
  forest_data, "sens", "sens_lo", "sens_hi", "pt_size_sens", "label_sens",
  xlab  = "Sensitivity (95% CI)",
  title = sprintf("Figure 2 - Forest Plot: Study-Level Sensitivities (rho_cluster = %.2f)",
                  primary_rho),
  color = "#1B4F8A",
  note  = sprintf("Range: %.3f-%.3f  |  No pooled estimate (threshold heterogeneity)",
                  min(forest_data$sens), max(forest_data$sens)))
print(fig2)
cat("Figure 2 saved.\n")

fig3 <- make_forest(
  forest_data, "spec", "spec_lo", "spec_hi", "pt_size_spec", "label_spec",
  xlab  = "Specificity (95% CI)",
  title = sprintf("Figure 3 - Forest Plot: Study-Level Specificities (rho_cluster = %.2f)",
                  primary_rho),
  color = "#2A9D8F",
  note  = sprintf("Range: %.3f-%.3f  |  No pooled estimate (threshold heterogeneity)",
                  min(forest_data$spec), max(forest_data$spec)))
print(fig3)
cat("Figure 3 saved.\n")

anchor_df <- data.frame(
  thresh = anchor_thresh_draws[valid],
  spec   = anchor_spec_draws[valid]
)

fig4a <- ggplot(anchor_df, aes(x = thresh * 100)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                 fill = "#2A9D8F", color = "white",
                 linewidth = 0.3, alpha = 0.8) +
  geom_density(color = "#1B4F8A", linewidth = 1) +
  geom_vline(xintercept = anchor_thresh_med * 100,
             color = "#2A9D8F", linewidth = 1.5) +
  geom_vline(xintercept = c(anchor_thresh_lo * 100, anchor_thresh_hi * 100),
             color = "#2A9D8F", linewidth = 1, linetype = "dashed") +
  annotate("text", x = anchor_thresh_med * 100 + 0.3, y = Inf, vjust = 1.5,
           label = sprintf("Median = %.2f%%\n95%% CrI: [%.2f%%, %.2f%%]",
                           anchor_thresh_med * 100, anchor_thresh_lo * 100,
                           anchor_thresh_hi * 100),
           color = "#2A9D8F", size = 3.5, hjust = 0) +
  labs(title    = sprintf("Posterior: dd-cfDNA Cutoff at %.0f%% Sensitivity (rho_cluster = %.2f)",
                          sens_target * 100, primary_rho),
       subtitle = sprintf("%d posterior draws | adaptive walk-down to %.0f%% target",
                          n_draws, sens_target * 100),
       x = sprintf("dd-cfDNA Threshold at %.0f%% Sensitivity (%%)",
                   sens_target * 100),
       y = "Posterior Density") +
  theme_bw(base_size = 11) +
  theme(panel.background = element_rect(fill = "#FAFAFA"),
        plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

fig4b <- ggplot(anchor_df, aes(x = spec)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                 fill = "#1B4F8A", color = "white",
                 linewidth = 0.3, alpha = 0.8) +
  geom_density(color = "#D94F3D", linewidth = 1) +
  geom_vline(xintercept = anchor_spec_med, color = "#1B4F8A", linewidth = 1.5) +
  geom_vline(xintercept = c(anchor_spec_lo, anchor_spec_hi),
             color = "#1B4F8A", linewidth = 1, linetype = "dashed") +
  annotate("text", x = anchor_spec_med + 0.02, y = Inf, vjust = 1.5,
           label = sprintf("Median = %.3f\n95%% CrI: [%.3f, %.3f]",
                           anchor_spec_med, anchor_spec_lo, anchor_spec_hi),
           color = "#1B4F8A", size = 3.5, hjust = 0) +
  labs(title    = sprintf("Posterior: Specificity at %.0f%% Sensitivity Anchor",
                          sens_target * 100),
       subtitle = "Derived from posterior threshold draws",
       x = sprintf("Specificity at %.0f%% Sensitivity Threshold",
                   sens_target * 100),
       y = "Posterior Density") +
  theme_bw(base_size = 11) +
  theme(panel.background = element_rect(fill = "#FAFAFA"),
        plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

fig4 <- fig4a / fig4b
print(fig4)
cat("Figure 4 saved.\n")

cat("\n=== Strategy-Based Decision Curve Analysis ===\n")

pt_range <- seq(0.01, 0.60, length.out = 300)
n_pt     <- length(pt_range)
n_strat_draws <- nrow(sens_mat)

bvn_cdf <- function(a, b, rho) {
  pbivnorm::pbivnorm(a, b, rep(rho, length(a)))
}

.copula_eps <- 1e-10
clip_marginal <- function(p) pmin(pmax(p, .copula_eps), 1 - .copula_eps)

copula_strat4_joint <- function(sens_vec, spec_vec, rho) {
  z_pos_dp <- qnorm(1 - clip_marginal(sens_vec))
  z_pos_dn <- qnorm(    clip_marginal(spec_vec))
  list(
    p_both_neg_dp = bvn_cdf(z_pos_dp, z_pos_dp, rho),
    p_both_neg_dn = bvn_cdf(z_pos_dn, z_pos_dn, rho)
  )
}

copula_strat5_joint <- function(sens_low_vec, sens_high_vec,
                                  spec_low_vec, spec_high_vec, rho) {
  n <- length(sens_low_vec)

  ok_order  <- (sens_low_vec  > sens_high_vec) &
                (spec_high_vec > spec_low_vec)
  ok_order[is.na(ok_order)] <- FALSE
  valid <- !is.na(sens_low_vec) & !is.na(sens_high_vec) &
           !is.na(spec_low_vec) & !is.na(spec_high_vec) &
           ok_order

  out <- list(
    p_high_dp     = sens_high_vec,
    p_low_dn      = spec_low_vec,
    p_bordhigh_dp = rep(NA_real_, n),
    p_bordbord_dp = rep(NA_real_, n),
    p_bordlow_dn  = rep(NA_real_, n),
    p_bordbord_dn = rep(NA_real_, n),
    valid         = valid
  )
  if (!any(valid)) return(out)

  s_low_v  <- clip_marginal(sens_low_vec[valid])
  s_high_v <- clip_marginal(sens_high_vec[valid])
  sp_low_v <- clip_marginal(spec_low_vec[valid])
  sp_high_v <- clip_marginal(spec_high_vec[valid])

  z_low_pos  <- qnorm(1 - s_low_v)
  z_high_pos <- qnorm(1 - s_high_v)
  F_LL_pos <- bvn_cdf(z_low_pos,  z_low_pos,  rho)
  F_LH_pos <- bvn_cdf(z_low_pos,  z_high_pos, rho)
  F_HH_pos <- bvn_cdf(z_high_pos, z_high_pos, rho)
  Phi_low_pos  <- pnorm(z_low_pos)
  Phi_high_pos <- pnorm(z_high_pos)
  out$p_bordhigh_dp[valid] <- (Phi_high_pos - F_HH_pos) -
                                 (Phi_low_pos  - F_LH_pos)
  out$p_bordbord_dp[valid] <- F_HH_pos - 2 * F_LH_pos + F_LL_pos

  z_low_neg  <- qnorm(sp_low_v)
  z_high_neg <- qnorm(sp_high_v)
  F_LL_neg <- bvn_cdf(z_low_neg,  z_low_neg,  rho)
  F_LH_neg <- bvn_cdf(z_low_neg,  z_high_neg, rho)
  F_HH_neg <- bvn_cdf(z_high_neg, z_high_neg, rho)
  out$p_bordlow_dn[valid]  <- F_LH_neg - F_LL_neg
  out$p_bordbord_dn[valid] <- F_HH_neg - 2 * F_LH_neg + F_LL_neg

  out
}

clip01 <- function(x) pmin(pmax(x, 0), 1)
assert_valid_prob <- function(x, label, rho, tol = 1e-6) {
  x_finite <- x[!is.na(x)]
  if (length(x_finite) == 0L) return(invisible())
  if (any(!is.finite(x_finite))) {
    stop(sprintf("Non-finite values in %s at rho_repeat = %.2f",
                 label, rho), call. = FALSE)
  }
  if (any(x_finite < -tol | x_finite > 1 + tol)) {
    bad <- range(x_finite)
    stop(sprintf(
      "Joint probability %s out of [0, 1] at rho_repeat = %.2f (range: [%.6e, %.6e])",
      label, rho, bad[1], bad[2]), call. = FALSE)
  }
  invisible()
}

strat4_eff <- function(j) {
  list(sens_eff = clip01(1 - j$p_both_neg_dp),
       spec_eff = clip01(j$p_both_neg_dn))
}

strat5_eff <- function(j, p_b) {
  list(
    sens_eff = clip01(j$p_high_dp + j$p_bordhigh_dp + p_b * j$p_bordbord_dp),
    spec_eff = clip01(j$p_low_dn  + j$p_bordlow_dn  + (1 - p_b) * j$p_bordbord_dn)
  )
}

net_benefit <- function(sens, spec, prev, pt) {
  sens * prev - (1 - spec) * (1 - prev) * (pt / (1 - pt))
}

prev_levels <- c(0.10, 0.20, 0.30)
cat(sprintf("rho_cluster (DEff) primary: %.2f\n", primary_rho))
cat(sprintf("rho_repeat  (DCA)  default: %.2f\n", rho_repeat_default))
cat("Borderline resolution: threshold-dependent p_b = 1 - p_t\n")
cat(sprintf("Prevalence levels: %s\n",
            paste(sprintf("%.0f%%", prev_levels * 100), collapse = ", ")))
cat("DCA joints: Gaussian copula on Bernoulli marginals (pbivnorm).\n")

j4 <- copula_strat4_joint(sens_youden_vec, spec_youden_vec,
                           rho_repeat_default)
assert_valid_prob(j4$p_both_neg_dp,
                   "p(both neg | D+, S4)", rho_repeat_default)
assert_valid_prob(j4$p_both_neg_dn,
                   "p(both neg | D-, S4)", rho_repeat_default)
eff4 <- strat4_eff(j4)
assert_valid_prob(eff4$sens_eff, "sens_eff (S4)", rho_repeat_default)
assert_valid_prob(eff4$spec_eff, "spec_eff (S4)", rho_repeat_default)

j5 <- copula_strat5_joint(sens_low_vec, sens_high_vec,
                           spec_low_vec, spec_high_vec,
                           rho_repeat_default)
valid_5 <- j5$valid
n_inverted_5 <- sum(!valid_5 & !is.na(sens_low_vec) & !is.na(sens_high_vec) &
                     !is.na(spec_low_vec) & !is.na(spec_high_vec))
cat(sprintf(
  "  Strategy 5 valid draws: %d / %d (%.1f%%) | order-inverted: %d | NA: %d\n",
  sum(valid_5), length(valid_5), 100 * mean(valid_5),
  n_inverted_5,
  sum(is.na(sens_low_vec) | is.na(sens_high_vec) |
      is.na(spec_low_vec) | is.na(spec_high_vec))))
for (nm in c("p_bordhigh_dp", "p_bordbord_dp",
              "p_bordlow_dn",  "p_bordbord_dn")) {
  assert_valid_prob(j5[[nm]][valid_5],
                     paste0(nm, " (S5)"), rho_repeat_default)
}

dca_all <- list()
for (prev in prev_levels) {
  prev_label <- sprintf("Prevalence = %.0f%%", prev * 100)
  cat(sprintf("  Computing strategies at %s ...\n", prev_label))

  nb_biopsy_all  <- matrix(NA, n_strat_draws, n_pt)
  nb_monitor_all <- matrix(NA, n_strat_draws, n_pt)
  nb_single      <- matrix(NA, n_strat_draws, n_pt)
  nb_repeat_neg  <- matrix(NA, n_strat_draws, n_pt)
  nb_borderline  <- matrix(NA, n_strat_draws, n_pt)

  for (p in seq_len(n_pt)) {
    nb_biopsy_all[, p]  <- net_benefit(1, 0, prev, pt_range[p])
    nb_monitor_all[, p] <- 0
    nb_single[, p]      <- net_benefit(sens_youden_vec, spec_youden_vec,
                                        prev, pt_range[p])
    nb_repeat_neg[, p]  <- net_benefit(eff4$sens_eff, eff4$spec_eff,
                                        prev, pt_range[p])

    p_b <- 1 - pt_range[p]
    s5  <- strat5_eff(j5, p_b)

    assert_valid_prob(s5$sens_eff[valid_5],
                       sprintf("sens_eff (S5, p_t = %.3f)", pt_range[p]),
                       rho_repeat_default)
    assert_valid_prob(s5$spec_eff[valid_5],
                       sprintf("spec_eff (S5, p_t = %.3f)", pt_range[p]),
                       rho_repeat_default)
    nb_borderline[valid_5, p] <- net_benefit(s5$sens_eff[valid_5],
                                              s5$spec_eff[valid_5],
                                              prev, pt_range[p])
  }

  summarise_nb <- function(mat) {
    list(med = apply(mat, 2, median, na.rm = TRUE),
         lo  = apply(mat, 2, quantile, 0.025, na.rm = TRUE),
         hi  = apply(mat, 2, quantile, 0.975, na.rm = TRUE))
  }
  s1 <- summarise_nb(nb_biopsy_all);  s2 <- summarise_nb(nb_monitor_all)
  s3 <- summarise_nb(nb_single);      s4 <- summarise_nb(nb_repeat_neg)
  s5 <- summarise_nb(nb_borderline)

  strategy_names <- c("Biopsy All", "Monitor All", "Single Test",
                       "Repeat-if-Negative", "Repeat-if-Borderline")
  ribbon_df <- data.frame(
    pt        = rep(pt_range, 5),
    median_nb = c(s1$med, s2$med, s3$med, s4$med, s5$med),
    lo        = c(s1$lo,  s2$lo,  s3$lo,  s4$lo,  s5$lo),
    hi        = c(s1$hi,  s2$hi,  s3$hi,  s4$hi,  s5$hi),
    strategy  = rep(strategy_names, each = n_pt),
    prev_label = prev_label
  )
  dca_all[[prev_label]] <- ribbon_df
}
dca_all_df <- do.call(rbind, dca_all)
dca_all_df$prev_label <- factor(dca_all_df$prev_label,
  levels = sprintf("Prevalence = %.0f%%", prev_levels * 100))
dca_all_df$strategy <- factor(dca_all_df$strategy,
  levels = c("Biopsy All", "Monitor All", "Single Test",
             "Repeat-if-Negative", "Repeat-if-Borderline"))

y_max <- max(prev_levels) + 0.03
ref_strategies  <- c("Biopsy All", "Monitor All")
test_strategies <- c("Single Test", "Repeat-if-Negative", "Repeat-if-Borderline")
dca_ref  <- dca_all_df[dca_all_df$strategy %in% ref_strategies, ]
dca_test <- dca_all_df[dca_all_df$strategy %in% test_strategies, ]

fig5 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3, alpha = 0.4) +
  geom_line(data = dca_ref,
            aes(x = pt, y = median_nb, color = strategy, linetype = strategy),
            linewidth = 0.8, na.rm = TRUE) +
  geom_ribbon(data = dca_test,
              aes(x = pt, ymin = lo, ymax = hi, fill = strategy),
              alpha = 0.12) +
  geom_line(data = dca_test,
            aes(x = pt, y = median_nb, color = strategy),
            linewidth = 1.1, na.rm = TRUE) +
  facet_wrap(~ prev_label, nrow = 1, scales = "free_y") +
  scale_color_manual(
    values = c("Biopsy All"           = "#DC3545",
               "Monitor All"          = "black",
               "Single Test"          = "#1B4F8A",
               "Repeat-if-Negative"   = "#E8850C",
               "Repeat-if-Borderline" = "#2A9D8F"),
    name = "Strategy") +
  scale_fill_manual(
    values = c("Single Test"          = "#1B4F8A",
               "Repeat-if-Negative"   = "#E8850C",
               "Repeat-if-Borderline" = "#2A9D8F"),
    name = "Strategy") +
  scale_linetype_manual(
    values = c("Biopsy All" = "longdash", "Monitor All" = "dotted"),
    name = "Strategy") +
  scale_x_continuous(labels = percent_format()) +
  coord_cartesian(ylim = c(-0.08, y_max)) +
  labs(title    = "Strategy-Based Decision Curve Analysis - dd-cfDNA",
       subtitle = sprintf(
         "5 strategies x 3 prevalences | rho_cluster = %.2f (primary), rho_repeat = %.2f | Posterior 95%% CrI",
         primary_rho, rho_repeat_default),
       x = "Threshold Probability (willingness to biopsy)",
       y = "Net Benefit",
       caption = sprintf(
         paste0("Single Test: dd-cfDNA at Youden threshold (%.2f%%) | ",
                "Repeat-if-Negative: retest in 3-7 days if negative | ",
                "Repeat-if-Borderline: three-zone (T_low=%.2f%%, T_high=%.2f%%) | ",
                "Repeat-test joints via Gaussian copula on Bernoulli marginals (pbivnorm)"),
         youden_thresh * 100, anchor_thresh_med * 100, spec90_thresh_med * 100)) +
  theme_bw(base_size = 11) +
  theme(panel.background = element_rect(fill = "#FAFAFA"),
        legend.position  = "bottom",
        legend.background = element_rect(color = "grey80", fill = "white"),
        legend.text   = element_text(size = 9),
        legend.key.width = unit(1.8, "cm"),
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        plot.caption  = element_text(size = 7.5, color = "grey50", hjust = 0),
        strip.text    = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "#EEF4FB", color = "grey80"))

print(fig5)
cat("Figure 5 (Strategy DCA) saved.\n")

cat("Computing rho_repeat sensitivity (rho = 0.30, 0.50, 0.70) ...\n")
rho_repeat_values <- c(0.30, 0.50, 0.70)
prev_rho <- 0.20
dca_rho_all <- list()

for (rho_r in rho_repeat_values) {
  rho_label <- sprintf("rho_repeat = %.2f", rho_r)

  j4_r <- copula_strat4_joint(sens_youden_vec, spec_youden_vec, rho_r)
  assert_valid_prob(j4_r$p_both_neg_dp, "p(both neg | D+, S4)", rho_r)
  assert_valid_prob(j4_r$p_both_neg_dn, "p(both neg | D-, S4)", rho_r)
  eff4_r <- strat4_eff(j4_r)
  assert_valid_prob(eff4_r$sens_eff, "sens_eff (S4)", rho_r)
  assert_valid_prob(eff4_r$spec_eff, "spec_eff (S4)", rho_r)

  j5_r <- copula_strat5_joint(sens_low_vec, sens_high_vec,
                                spec_low_vec, spec_high_vec, rho_r)

  valid_5_r <- j5_r$valid
  for (nm in c("p_bordhigh_dp", "p_bordbord_dp",
                "p_bordlow_dn",  "p_bordbord_dn")) {
    assert_valid_prob(j5_r[[nm]][valid_5_r], paste0(nm, " (S5)"), rho_r)
  }

  nb_single_rho     <- matrix(NA, n_strat_draws, n_pt)
  nb_repeat_neg_rho <- matrix(NA, n_strat_draws, n_pt)
  nb_borderline_rho <- matrix(NA, n_strat_draws, n_pt)
  for (p in seq_len(n_pt)) {
    nb_single_rho[, p] <- net_benefit(sens_youden_vec, spec_youden_vec,
                                       prev_rho, pt_range[p])
    nb_repeat_neg_rho[, p] <- net_benefit(eff4_r$sens_eff, eff4_r$spec_eff,
                                            prev_rho, pt_range[p])
    p_b <- 1 - pt_range[p]
    s5  <- strat5_eff(j5_r, p_b)
    assert_valid_prob(s5$sens_eff[valid_5_r],
                       sprintf("sens_eff (S5, p_t = %.3f)", pt_range[p]),
                       rho_r)
    assert_valid_prob(s5$spec_eff[valid_5_r],
                       sprintf("spec_eff (S5, p_t = %.3f)", pt_range[p]),
                       rho_r)
    nb_borderline_rho[valid_5_r, p] <- net_benefit(s5$sens_eff[valid_5_r],
                                                    s5$spec_eff[valid_5_r],
                                                    prev_rho, pt_range[p])
  }
  med_only <- function(mat) apply(mat, 2, median, na.rm = TRUE)
  df_rho <- data.frame(
    pt        = rep(pt_range, 3),
    median_nb = c(med_only(nb_single_rho),
                  med_only(nb_repeat_neg_rho),
                  med_only(nb_borderline_rho)),
    strategy  = rep(c("Single Test", "Repeat-if-Negative", "Repeat-if-Borderline"),
                    each = n_pt),
    rho_label = rho_label
  )
  dca_rho_all[[rho_label]] <- df_rho
}
dca_rho_df <- do.call(rbind, dca_rho_all)
dca_rho_df$rho_label <- factor(dca_rho_df$rho_label,
  levels = sprintf("rho_repeat = %.2f", rho_repeat_values))
dca_rho_df$strategy <- factor(dca_rho_df$strategy,
  levels = c("Single Test", "Repeat-if-Negative", "Repeat-if-Borderline"))

fig5b <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3, alpha = 0.4) +
  geom_line(data = dca_rho_df,
            aes(x = pt, y = median_nb, color = strategy),
            linewidth = 1.0, na.rm = TRUE) +
  facet_wrap(~ rho_label, nrow = 1) +
  scale_color_manual(
    values = c("Single Test"          = "#1B4F8A",
               "Repeat-if-Negative"   = "#E8850C",
               "Repeat-if-Borderline" = "#2A9D8F"),
    name = "Strategy") +
  scale_x_continuous(labels = percent_format()) +
  coord_cartesian(ylim = c(-0.05, prev_rho + 0.03)) +
  labs(title    = "rho_repeat Sensitivity - Repeat-Testing Strategies",
       subtitle = sprintf("Prevalence = %.0f%% | varying intra-patient repeat correlation",
                          prev_rho * 100),
       x = "Threshold Probability", y = "Net Benefit",
       caption = "Higher rho_repeat = more correlated repeat tests (less independent information)") +
  theme_bw(base_size = 11) +
  theme(panel.background = element_rect(fill = "#FAFAFA"),
        legend.position  = "bottom",
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        plot.caption  = element_text(size = 8, color = "grey50", hjust = 0),
        strip.text    = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "#EEF4FB", color = "grey80"))

print(fig5b)
cat("Figure 5b saved.\n")

T1 <- anchor_thresh_med
T2 <- spec90_thresh_med
T_youden <- youden_thresh
x_max_pct <- ceiling(max(df$threshold) * 100 * 1.25)

zone_df <- data.frame(
  xmin = c(0, T1, T2) * 100,
  xmax = c(T1, T2, x_max_pct),
  zone = c("Rule-Out", "Indeterminate", "Rule-In")
)

fig6 <- ggplot() +
  geom_rect(data = zone_df,
            aes(xmin = xmin, xmax = xmax, ymin = 0.20, ymax = 0.90,
                fill = zone, color = zone),
            linewidth = 1, alpha = 0.80) +
  scale_fill_manual(values  = c("Rule-Out"      = "#D4EDDA",
                                "Indeterminate" = "#FFF3CD",
                                "Rule-In"       = "#F8D7DA")) +
  scale_color_manual(values = c("Rule-Out"      = "#28A745",
                                "Indeterminate" = "#FFC107",
                                "Rule-In"       = "#DC3545")) +
  annotate("text", x = T1 * 100 / 2, y = 0.84,
           label = "ZONE 1: RULE-OUT",
           fontface = "bold", size = 3.8, color = "#155724") +
  annotate("text", x = (T1 + T2) * 100 / 2, y = 0.84,
           label = "ZONE 2:\nINDETERMINATE",
           fontface = "bold", size = 3.2, color = "#856404") +
  annotate("text", x = (T2 * 100 + x_max_pct) / 2, y = 0.84,
           label = "ZONE 3: RULE-IN",
           fontface = "bold", size = 3.8, color = "#721C24") +
  annotate("text", x = T1 * 100 / 2, y = 0.76,
           label = sprintf("dd-cfDNA < %.1f%%", T1 * 100),
           size = 3, color = "#155724", fontface = "italic") +
  annotate("text", x = (T1 + T2) * 100 / 2, y = 0.76,
           label = sprintf("%.1f%% <= dd-cfDNA\n<= %.1f%%", T1 * 100, T2 * 100),
           size = 2.8, color = "#856404", fontface = "italic") +
  annotate("text", x = (T2 * 100 + x_max_pct) / 2, y = 0.76,
           label = sprintf("dd-cfDNA > %.1f%%", T2 * 100),
           size = 3, color = "#721C24", fontface = "italic") +
  annotate("text", x = T1 * 100 / 2, y = 0.62,
           label = sprintf("Sensitivity ~ %.0f%%\nFNR < %.0f%% | High NPV",
                           sens_target * 100, (1 - sens_target) * 100),
           size = 3, color = "#155724", lineheight = 1.5) +
  annotate("text", x = (T1 + T2) * 100 / 2, y = 0.62,
           label = "Uncertain Sens/Spec\nIntermediate zone",
           size = 3, color = "#856404", lineheight = 1.5) +
  annotate("text", x = (T2 * 100 + x_max_pct) / 2, y = 0.62,
           label = sprintf("Spec >= %.0f%% | FPR <= %.0f%%\nSens = %.0f%% [%.0f-%.0f%%]",
                           spec_target * 100, (1 - spec_target) * 100,
                           spec90_sens_med * 100, spec90_sens_lo * 100,
                           spec90_sens_hi * 100),
           size = 3, color = "#721C24", lineheight = 1.5) +
  geom_label(data = data.frame(
               x = T1 * 100 / 2, y = 0.37,
               lab = "Rejection unlikely\n-> Avoid biopsy\n-> Routine surveillance"),
             aes(x = x, y = y, label = lab),
             size = 2.8, color = "#155724", fill = "#D4EDDA",
             label.size = 0.4, lineheight = 1.4, inherit.aes = FALSE) +
  geom_label(data = data.frame(
               x = (T1 + T2) * 100 / 2, y = 0.37,
               lab = "Uncertain -> Combine with:\n- Clinical symptoms\n- Echo / GEP / biomarkers"),
             aes(x = x, y = y, label = lab),
             size = 2.6, color = "#856404", fill = "#FFF3CD",
             label.size = 0.4, lineheight = 1.4, inherit.aes = FALSE) +
  geom_label(data = data.frame(
               x = (T2 * 100 + x_max_pct) / 2, y = 0.37,
               lab = "High suspicion for AR\n-> Proceed to EMB\n-> Intensify surveillance"),
             aes(x = x, y = y, label = lab),
             size = 2.8, color = "#721C24", fill = "#F8D7DA",
             label.size = 0.4, lineheight = 1.4, inherit.aes = FALSE) +
  geom_vline(xintercept = T1 * 100, color = "#2A9D8F",
             linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = T2 * 100, color = "#7B2D8E",
             linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = T_youden * 100, color = "#E8850C",
             linewidth = 0.8, linetype = "dotted") +
  annotate("text", x = T1 * 100, y = 0.94,
           label = sprintf("T1 = %.1f%%\n(%s)", T1 * 100, anchor_sens_label),
           size = 3, color = "#2A9D8F", fontface = "bold") +
  annotate("text", x = T2 * 100, y = 0.94,
           label = sprintf("T2 = %.1f%%\n(%s)", T2 * 100, anchor_spec_label),
           size = 3, color = "#7B2D8E", fontface = "bold") +
  annotate("text", x = T_youden * 100, y = 0.17,
           label = sprintf("Youden = %.1f%%\n(J = %.2f)", T_youden * 100, youden_J),
           size = 2.5, color = "#E8850C", fontface = "italic") +
  annotate("text", x = x_max_pct / 2, y = 0.10,
           label = sprintf(
             "T1 = %.2f%% (CrI: %.2f-%.2f%%) | Spec at T1 = %.3f [%.3f-%.3f]   |   T2 = %.2f%% (CrI: %.2f-%.2f%%) | Sens at T2 = %.3f [%.3f-%.3f]   |   Youden = %.2f%%   |   rho_cluster = %.2f\nSummary point (marginalised over study thresholds): Sens = %.3f [%.3f-%.3f], Spec = %.3f [%.3f-%.3f]   |   Operating point at mean observed threshold (%s%.3f): Sens = %.3f, Spec = %.3f",
             T1 * 100, anchor_thresh_lo * 100, anchor_thresh_hi * 100,
             anchor_spec_med, anchor_spec_lo, anchor_spec_hi,
             T2 * 100, spec90_thresh_lo * 100, spec90_thresh_hi * 100,
             spec90_sens_med, spec90_sens_lo, spec90_sens_hi,
             T_youden * 100, primary_rho,
             summary_sens, summary_sens_lo, summary_sens_hi,
             summary_spec, summary_spec_lo, summary_spec_hi,
             "γ₀ = ", median(gamma_0_draws),
             gamma0_sens, gamma0_spec),
           size = 2.0, color = "grey50", fontface = "italic") +
  scale_x_continuous(limits = c(0, x_max_pct),
                     breaks = unique(round(c(0, T1 * 100, T_youden * 100,
                                              T2 * 100,
                                              seq(0, x_max_pct, length.out = 5)),
                                            1)),
                     labels = function(x) sprintf("%.1f%%", x)) +
  scale_y_continuous(limits = c(0.05, 1.02)) +
  labs(title    = "Clinical Decision Map - dd-cfDNA Threshold Framework",
       subtitle = sprintf(
         "T1 = %s (rule-out) | T2 = %s (rule-in) | rho_cluster = %.2f (primary)",
         anchor_sens_label, anchor_spec_label, primary_rho),
       x = "dd-cfDNA (%)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.background = element_blank(),
        panel.grid       = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text.y      = element_blank(),
        axis.text.x      = element_text(size = 8),
        legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 11),
        plot.subtitle    = element_text(size = 9, color = "grey40"))

print(fig6)
cat("Figure 6 saved.\n")

cat("\n============================================================\n")
cat("FINAL RESULTS SUMMARY\n")
cat("============================================================\n")
cat(sprintf("Primary rho_cluster:              %.2f\n", primary_rho))
cat(sprintf("rho sweep:                        {%s}\n",
            paste(sprintf("%.2f", rho_values), collapse = ", ")))
cat(sprintf("Mean DEff at primary:             %.2f\n",
            mean(primary$df_eff$deff)))
cat(sprintf("Total effective N at primary:     %d\n",
            sum(primary$df_eff$n_pos_eff + primary$df_eff$n_neg_eff)))

cat("\n--- Primary results (rho = ", sprintf("%.2f", primary_rho),
    ") ---\n", sep = "")
cat(sprintf("HSROC AUC:                        %.3f (95%% CrI: %.3f-%.3f)\n",
            auc_sroc, auc_lo, auc_hi))
cat("Summary point (marginalised over empirical threshold distribution):\n")
cat(sprintf("  Sensitivity:                    %.3f (95%% CrI: %.3f-%.3f)\n",
            summary_sens, summary_sens_lo, summary_sens_hi))
cat(sprintf("  Specificity:                    %.3f (95%% CrI: %.3f-%.3f)\n",
            summary_spec, summary_spec_lo, summary_spec_hi))
cat("Operating point at the mean observed threshold (gamma_0):\n")
cat(sprintf("  Sensitivity:                    %.3f (95%% CrI: %.3f-%.3f)\n",
            gamma0_sens, gamma0_sens_lo, gamma0_sens_hi))
cat(sprintf("  Specificity:                    %.3f (95%% CrI: %.3f-%.3f)\n",
            gamma0_spec, gamma0_spec_lo, gamma0_spec_hi))
cat(sprintf("Youden cutoff:                    %.3f%%\n", youden_thresh * 100))
cat(sprintf("Youden sensitivity:               %.3f\n", youden_sens))
cat(sprintf("Youden specificity:               %.3f\n", youden_spec))
cat(sprintf("Youden J:                         %.3f [%.3f-%.3f]\n",
            youden_J, youden_J_lo, youden_J_hi))
cat(sprintf("%.0f%% Sens anchor cutoff:        %.3f%% (CrI: %.3f%%-%.3f%%)\n",
            sens_target * 100, anchor_thresh_med * 100,
            anchor_thresh_lo * 100, anchor_thresh_hi * 100))
cat(sprintf("Specificity at sens anchor:       %.3f (CrI: %.3f-%.3f)\n",
            anchor_spec_med, anchor_spec_lo, anchor_spec_hi))
cat(sprintf("%.0f%% Spec anchor cutoff:        %.3f%% (CrI: %.3f%%-%.3f%%)\n",
            spec_target * 100, spec90_thresh_med * 100,
            spec90_thresh_lo * 100, spec90_thresh_hi * 100))
cat(sprintf("Sensitivity at spec anchor:       %.3f (CrI: %.3f-%.3f)\n",
            spec90_sens_med, spec90_sens_lo, spec90_sens_hi))
cat(sprintf("PSIS+exact blended elpd_loo:      %.2f (se = %.2f)\n",
            elpd_loo_total, elpd_loo_total_se))
cat(sprintf("p_loo (blended; in-sample - elpd_loo): %.2f | p_loo/N = %.3f\n",
            p_loo_blended, p_loo_blended / N))
cat(sprintf("p_loo (pure PSIS, from loo::loo):   %.2f (se = %.2f) | p_loo/N = %.3f\n",
            p_loo_total, p_loo_total_se, p_loo_total / N))
cat(sprintf("LOOIC:                            %.2f (se = %.2f)\n",
            looic, looic_se))
cat(sprintf("max Pareto-k | studies k > 0.7:   %.2f | %d / %d\n",
            pareto_k_max, n_high_k, N))
if (n_high_k > 0L) {
  cat(sprintf("  Exact-refit substitutions made for: %s\n",
              paste(loo_study[!is.na(exact_lpd_per_study)], collapse = ", ")))
}

cat("\n--- Sweep ranges (across rho in {0, 0.30, 0.50, 0.70}) ---\n")
fmt_range <- function(vals, fmt = "%.3f") {
  sprintf("min %s | median %s | max %s",
          sprintf(fmt, min(vals)),
          sprintf(fmt, median(vals)),
          sprintf(fmt, max(vals)))
}
auc_meds  <- sapply(results, function(r) r$auc_sroc)
yth_pcts  <- sapply(results, function(r) r$youden_thresh * 100)
yj_meds   <- sapply(results, function(r) r$youden_J)
sens_tgts <- sapply(results, function(r) r$sens_target)
spec_tgts <- sapply(results, function(r) r$spec_target)
sa_pcts   <- sapply(results, function(r) r$anchor_thresh_med * 100)
spa_pcts  <- sapply(results, function(r) r$spec90_thresh_med * 100)
elpd_vals <- sapply(results, function(r) r$elpd_loo_total)
ploo_vals <- sapply(results, function(r) r$p_loo_total)
ploo_blended_vals <- sapply(results, function(r) r$p_loo_blended)
kmax_vals <- sapply(results, function(r) r$pareto_k_max)
nhk_vals  <- sapply(results, function(r) r$n_high_k)

cat(sprintf("AUC across rho:                   %s\n", fmt_range(auc_meds)))
cat(sprintf("Youden cutoff (%%) across rho:     %s\n",
            fmt_range(yth_pcts, "%.3f%%")))
cat(sprintf("Youden J across rho:              %s\n", fmt_range(yj_meds)))
cat(sprintf("Sens-anchor target (%%) across rho: %s\n",
            paste(sprintf("%.0f%%", sens_tgts * 100), collapse = " | ")))
cat(sprintf("Sens-anchor cutoff (%%) across rho: %s\n",
            fmt_range(sa_pcts, "%.3f%%")))
cat(sprintf("Spec-anchor target (%%) across rho: %s\n",
            paste(sprintf("%.0f%%", spec_tgts * 100), collapse = " | ")))
cat(sprintf("Spec-anchor cutoff (%%) across rho: %s\n",
            fmt_range(spa_pcts, "%.3f%%")))
cat(sprintf("elpd_loo across rho:              %s\n", fmt_range(elpd_vals, "%.2f")))
cat(sprintf("p_loo (blended) across rho:       %s\n", fmt_range(ploo_blended_vals, "%.2f")))
cat(sprintf("p_loo (pure PSIS) across rho:     %s\n", fmt_range(ploo_vals, "%.2f")))
cat(sprintf("max Pareto-k across rho:          %s\n", fmt_range(kmax_vals, "%.2f")))
cat(sprintf("# studies k > 0.7 across rho:     %s\n",
            paste(nhk_vals, collapse = " | ")))
