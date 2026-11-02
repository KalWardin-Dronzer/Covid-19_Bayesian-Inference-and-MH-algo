# Bayesian SIR fit to India's first COVID-19 cases (14-31 March 2020) by adaptive Metropolis-Hastings.
#
#   Rscript covid_sir.R        -> results/summary.csv, results/RESULTS.md, results/fit.png
#
# What 18 days of early growth can identify is the growth rate, r = beta - gamma; beta and gamma separately, and so
# R0 = beta / gamma, depend on what is assumed about the infectious period. The script therefore reports r and the
# doubling time as the estimates, and R0 under three priors on gamma.

suppressPackageStartupMessages({
  library(deSolve)
  library(coda)
  library(ggplot2)
  library(patchwork)
})

# -- Data -------------------------------------------------------------------------------------------------------
# Johns Hopkins CSSE COVID-19 time series, confirmed cases, India (CSSEGISandData/COVID-19,
# csse_covid_19_time_series/time_series_covid19_confirmed_global.csv). New cases are daily differences.
data <- read.csv("data/india_confirmed_jhu_march2020.csv")
obs <- data$new_cases
dates <- as.Date(data$date)
T_len <- length(obs)
N <- 1.38e9

# -- Model ------------------------------------------------------------------------------------------------------
sir <- function(t, y, p) {
  with(as.list(c(y, p)), {
    infections <- beta * S * I / N
    list(c(-infections, infections - gamma * I))
  })
}

incidence <- function(beta, gamma, I0) {
  out <- ode(c(S = N - I0, I = I0), 0:T_len, sir, c(beta = beta, gamma = gamma), method = "lsoda")
  pmax(-diff(out[, "S"]), 1e-9)
}

# Parameters are sampled as theta = (r, log gamma, log I0, log k), where r = beta - gamma is the early growth rate
# and k the negative-binomial dispersion. In (beta, gamma) the posterior is a narrow ridge along beta - gamma =
# r, which the original random walk could barely move along (effective sample size about 19 in 30,000 draws).
unpack <- function(theta) {
  gamma <- exp(theta[2])
  c(beta = theta[1] + gamma, gamma = gamma, I0 = exp(theta[3]), k = exp(theta[4]))
}

log_post <- function(theta, prior) {
  p <- unpack(theta)
  if (p["beta"] <= 0 || p["gamma"] <= 0) return(-Inf)
  R0 <- p["beta"] / p["gamma"]
  if (R0 < 1 || R0 > 8) return(-Inf)
  lp <- dnorm(p["beta"], 0.35, 0.08, log = TRUE) +                  # beta prior, as in the original
    dnorm(p["gamma"], prior$mean, prior$sd, log = TRUE) +           # gamma prior: varied below
    theta[2] +                                                      # Jacobian of (r, log gamma) -> (beta, gamma)
    dnorm(theta[3], log(50), 1.5, log = TRUE) +                     # infections on 13 March: weakly informative
    dnorm(theta[4], log(10), 1.5, log = TRUE)                       # dispersion: weakly informative
  mu <- incidence(p["beta"], p["gamma"], p["I0"])
  unname(lp + sum(dnbinom(obs, size = p["k"], mu = mu, log = TRUE)))
}

# -- Adaptive Metropolis-Hastings ---------------------------------------------------------------------------------
# Gaussian random walk; during burn-in the proposal covariance is re-estimated every 500 iterations from the chain
# so far (scaled by 2.38^2 / d, Haario et al. 2001), then frozen, so the kept draws come from a fixed kernel.
run_chain <- function(start, prior, n_iter, burn_in) {
  d <- length(start)
  cov_prop <- diag(c(0.01, 0.05, 0.1, 0.1)^2)
  draws <- matrix(NA_real_, n_iter, d)
  current <- start
  current_lp <- log_post(current, prior)
  accepted <- 0
  for (i in 1:n_iter) {
    proposal <- current + drop(rnorm(d) %*% chol(cov_prop))
    proposal_lp <- log_post(proposal, prior)
    if (log(runif(1)) < proposal_lp - current_lp) {
      current <- proposal
      current_lp <- proposal_lp
      if (i > burn_in) accepted <- accepted + 1
    }
    draws[i, ] <- current
    if (i < burn_in && i >= 1000 && i %% 500 == 0) {
      cov_prop <- (2.38^2 / d) * cov(draws[(i %/% 2):i, ]) + diag(1e-8, d)
    }
  }
  kept <- draws[(burn_in + 1):n_iter, ]
  colnames(kept) <- c("r", "log_gamma", "log_I0", "log_k")
  list(draws = kept, acceptance = accepted / (n_iter - burn_in))
}

fit <- function(prior, chains = 4, n_iter = 20000, burn_in = 8000, seed = 2020) {
  set.seed(seed)
  start_values <- lapply(1:chains, function(j) c(runif(1, 0.10, 0.30), log(prior$mean) + rnorm(1, 0, 0.2),
                                                 log(50) + rnorm(1, 0, 0.5), log(10) + rnorm(1, 0, 0.5)))
  runs <- lapply(start_values, run_chain, prior = prior, n_iter = n_iter, burn_in = burn_in)
  list(chains = mcmc.list(lapply(runs, function(x) mcmc(x$draws))),
       acceptance = sapply(runs, `[[`, "acceptance"))
}

derived <- function(chains) {
  m <- as.matrix(chains)
  gamma <- exp(m[, "log_gamma"])
  beta <- m[, "r"] + gamma
  data.frame(r = m[, "r"], doubling_days = log(2) / m[, "r"], beta = beta, gamma = gamma,
             infectious_days = 1 / gamma, R0 = beta / gamma, I0 = exp(m[, "log_I0"]), k = exp(m[, "log_k"]))
}

ci <- function(x) sprintf("%.2f [%.2f, %.2f]", mean(x), quantile(x, 0.025), quantile(x, 0.975))

# -- Fits -------------------------------------------------------------------------------------------------------
priors <- list(
  original = list(mean = 0.10, sd = 0.025, label = "N(0.10, 0.025²), the original (10-day infectious period)"),
  shorter = list(mean = 0.20, sd = 0.05, label = "N(0.20, 0.05²) (5 days)"),
  longer = list(mean = 0.07, sd = 0.02, label = "N(0.07, 0.02²) (14 days)")
)
dir.create("results", showWarnings = FALSE)
summary_rows <- list()
fits <- list()
for (name in names(priors)) {
  started <- Sys.time()
  f <- fit(priors[[name]])
  fits[[name]] <- f
  d <- derived(f$chains)
  rhat <- gelman.diag(f$chains, autoburnin = FALSE)$psrf[, 1]
  ess <- effectiveSize(f$chains)
  summary_rows[[name]] <- data.frame(
    prior = name, r_mean = mean(d$r), r_lo = quantile(d$r, 0.025), r_hi = quantile(d$r, 0.975),
    doubling_mean = mean(d$doubling_days), doubling_lo = quantile(d$doubling_days, 0.025),
    doubling_hi = quantile(d$doubling_days, 0.975), R0_mean = mean(d$R0), R0_lo = quantile(d$R0, 0.025),
    R0_hi = quantile(d$R0, 0.975), infectious_mean = mean(d$infectious_days), I0_mean = mean(d$I0),
    k_mean = mean(d$k), cor_beta_gamma = cor(d$beta, d$gamma), max_rhat = max(rhat), min_ess = min(ess),
    acceptance = mean(f$acceptance), minutes = as.numeric(difftime(Sys.time(), started, units = "mins")))
  cat(sprintf("%-9s r %s/day, doubling %s days, R0 %s, max R-hat %.3f, min ESS %.0f\n", name, ci(d$r),
              ci(d$doubling_days), ci(d$R0), max(rhat), min(ess)))
}
summary <- do.call(rbind, summary_rows)
write.csv(summary, "results/summary.csv", row.names = FALSE)

# -- Posterior predictive check (original prior) --------------------------------------------------------------------
d <- derived(fits$original$chains)
set.seed(7)
pick <- sample(nrow(d), 500)
sims <- t(sapply(pick, function(i) rnbinom(T_len, size = d$k[i], mu = incidence(d$beta[i], d$gamma[i], d$I0[i]))))
band <- data.frame(date = dates, observed = obs, mean = colMeans(sims),
                   lo = apply(sims, 2, quantile, 0.025), hi = apply(sims, 2, quantile, 0.975))
coverage <- mean(obs >= band$lo & obs <= band$hi)

p1 <- ggplot(band, aes(date)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#9fc2ec") +
  geom_line(aes(y = mean), colour = "#2a78d6") +
  geom_point(aes(y = observed), colour = "#c0392b") +
  geom_vline(xintercept = as.Date("2020-03-25"), linetype = "dashed", colour = "grey40") +
  labs(title = "Posterior predictive: new cases per day (95% band)", x = NULL, y = "new confirmed cases") +
  theme_minimal(base_size = 10)
ridge <- do.call(rbind, lapply(names(fits), function(n) {
  x <- derived(fits[[n]]$chains)
  x[sample(nrow(x), 1500), c("beta", "gamma")] |> transform(prior = n)
}))
p2 <- ggplot(ridge, aes(gamma, beta, colour = prior)) +
  geom_point(size = 0.4, alpha = 0.5) +
  labs(title = "Posterior draws under three priors on gamma", x = "gamma (recovery rate)", y = "beta") +
  theme_minimal(base_size = 10)
ggsave("results/fit.png", p1 + p2, width = 11, height = 4, dpi = 130)

# -- Report ---------------------------------------------------------------------------------------------------
fmt <- function(m, lo, hi, digits = 2) sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), m, lo, hi)
rows <- sapply(names(priors), function(n) {
  s <- summary[summary$prior == n, ]
  sprintf("| %s | %s | %s | %s | %.1f | %.2f | %.3f | %.0f |", priors[[n]]$label,
          fmt(s$r_mean, s$r_lo, s$r_hi, 3), fmt(s$doubling_mean, s$doubling_lo, s$doubling_hi, 1),
          fmt(s$R0_mean, s$R0_lo, s$R0_hi), s$infectious_mean, s$cor_beta_gamma, s$max_rhat, s$min_ess)
})
report <- c(
  "# Results", "",
  "Generated by `Rscript covid_sir.R`. Do not edit by hand.", "",
  sprintf("Daily new confirmed cases in India, %s to %s (Johns Hopkins CSSE), SIR model with a negative-binomial likelihood, the number infectious on 13 March estimated. 4 chains of 20,000 adaptive Metropolis-Hastings iterations each, the first 8,000 discarded.", min(dates), max(dates)),
  "",
  "| Prior on gamma | Growth rate r per day | Doubling time (days) | R0 | Infectious period (days) | corr(beta, gamma) | Max R-hat | Min ESS |",
  "|---|---|---|---|---|---|---|---|",
  rows, "",
  "Posterior means with 95% credible intervals. The growth rate and doubling time barely move between priors; R0 follows the assumed infectious period, because 18 days of early growth identify beta - gamma but not beta and gamma separately.",
  "",
  sprintf("Posterior predictive check (original prior, 500 draws): %d of %d observed days fall inside the 95%% band.", sum(obs >= band$lo & obs <= band$hi), T_len),
  sprintf("Acceptance rate after burn-in: %.0f%% (original prior).", 100 * summary$acceptance[summary$prior == "original"]),
  sprintf("R %s, deSolve %s, coda %s.", getRversion(), packageVersion("deSolve"), packageVersion("coda")), "")
writeLines(report, "results/RESULTS.md")
cat("wrote results/RESULTS.md, results/summary.csv and results/fit.png\n")
