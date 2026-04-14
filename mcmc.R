# MCMC for India COVID-19 — Pre-lockdown phase (March 14-31, 2020)
# 18 days of unmitigated spread. Shows SIR parameter identifiability.

library(deSolve)
library(coda)
library(ggplot2)
library(patchwork)

# ── Real India COVID-19 data: March 14–31, 2020 (pre-lockdown / early lockdown) ─
# India lockdown announced March 24, nationwide from March 25
obs <- c(3, 5, 10, 13, 18, 24, 32, 42, 56, 73, 83, 87, 90, 109, 117, 148, 176, 197)
T_len <- length(obs)
N <- 1380000000
I0 <- 3  # seed: 3 cases on March 14

# ── SIR model ──────────────────────────────────────────────────────────────────
sir_ode <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

simulate <- function(beta, gamma) {
  init <- c(S = N - I0, I = I0, R = 0.0)
  times <- seq(0, T_len, by = 1)
  params <- c(beta = beta, gamma = gamma)
  
  out <- tryCatch(
    ode(y = init, times = times, func = sir_ode, parms = params, method = "lsoda"),
    error = function(e) NULL
  )
  
  if (is.null(out)) return(NULL)
  
  diff_S <- -diff(out[, "S"])
  return(pmax(diff_S, 0)) # new_cases = -ΔS
}

# ── Log-posterior ──────────────────────────────────────────────────────────────
sigma <- 8.0  # observation noise (counts are small here)

log_lik <- function(pred, obs) {
  if (is.null(pred) || any(pred <= 0)) return(-Inf)
  return(-0.5 * sum(((obs - pred) / sigma)^2))
}

log_prior <- function(beta, gamma) {
  # Informative priors based on COVID-19 epidemiology
  # beta ~ N(0.35, 0.08^2),  gamma ~ N(0.10, 0.025^2)
  if (beta <= 0 || gamma <= 0) return(-Inf)
  
  R0 <- beta / gamma
  if (R0 < 1.0 || R0 > 8.0) return(-Inf)
  
  lp_b <- -0.5 * ((beta - 0.35) / 0.08)^2
  lp_g <- -0.5 * ((gamma - 0.10) / 0.025)^2
  return(lp_b + lp_g)
}

log_post <- function(beta, gamma) {
  lp <- log_prior(beta, gamma)
  if (!is.finite(lp)) return(-Inf)
  pred <- simulate(beta, gamma)
  return(lp + log_lik(pred, obs))
}

# ── MAP via grid ───────────────────────────────────────────────────────────────
cat("Step 1: Grid search for MAP...\n")
best_lp <- -Inf
best_b <- 0.35
best_g <- 0.10

beta_grid <- seq(0.15, 0.65, length.out = 25)
gamma_grid <- seq(0.04, 0.25, length.out = 25)

for (b in beta_grid) {
  for (g in gamma_grid) {
    lp <- log_post(b, g)
    if (lp > best_lp) {
      best_lp <- lp
      best_b <- b
      best_g <- g
    }
  }
}
cat(sprintf("   MAP: β=%.4f, γ=%.4f, R₀=%.2f\n", best_b, best_g, best_b/best_g))

# ── Metropolis-Hastings ────────────────────────────────────────────────────────
cat("\nStep 2: Metropolis-Hastings (40,000 iterations)...\n")
set.seed(7)

n_iter <- 40000
burn_in <- 10000
prop_sd <- c(0.010, 0.004)

beta_c <- best_b
gamma_c <- best_g
lp_c <- best_lp

chain_b <- numeric(n_iter)
chain_g <- numeric(n_iter)
accepted <- 0

for (i in 1:n_iter) {
  bp <- beta_c + rnorm(1, 0, prop_sd[1])
  gp <- gamma_c + rnorm(1, 0, prop_sd[2])
  lp_p <- log_post(bp, gp)
  
  if (log(runif(1) + 1e-300) < (lp_p - lp_c)) {
    beta_c <- bp
    gamma_c <- gp
    lp_c <- lp_p
    accepted <- accepted + 1
  }
  
  chain_b[i] <- beta_c
  chain_g[i] <- gamma_c
  
  # Adapt proposal during burn-in
  if (i > 500 && i < burn_in && i %% 500 == 0) {
    ar <- accepted / i
    prop_sd <- prop_sd * exp(0.2 * (ar - 0.234))
    prop_sd <- pmax(pmin(prop_sd, 0.2), 1e-4)
  }
  
  if (i %% 10000 == 0) {
    ar <- accepted / i
    cat(sprintf("   Iter %5d | β=%.5f γ=%.5f R₀=%.3f | Accept: %.1f%%\n", 
                i, beta_c, gamma_c, beta_c/gamma_c, ar * 100))
  }
}

# Post-burn-in processing
post_b <- chain_b[(burn_in + 1):n_iter]
post_g <- chain_g[(burn_in + 1):n_iter]
post_R0 <- post_b / post_g

cat(sprintf("\n✅ Acceptance rate: %.1f%%\n", (accepted/n_iter) * 100))
ci_b <- quantile(post_b, c(0.025, 0.975))
cat(sprintf("β  = %.5f  [95%% CI: %.5f–%.5f]\n", mean(post_b), ci_b[1], ci_b[2]))

ci_g <- quantile(post_g, c(0.025, 0.975))
cat(sprintf("γ  = %.5f  [95%% CI: %.5f–%.5f]\n", mean(post_g), ci_g[1], ci_g[2]))

ci_R0 <- quantile(post_R0, c(0.025, 0.975))
cat(sprintf("R₀ = %.3f    [95%% CI: %.3f–%.3f]\n", mean(post_R0), ci_R0[1], ci_R0[2]))

cat(sprintf("ESS(β)=%.0f  ESS(γ)=%.0f\n", effectiveSize(post_b), effectiveSize(post_g)))

# ── Posterior predictive ───────────────────────────────────────────────────────
idx <- sample(length(post_b), 500, replace = FALSE)
curves <- list()
for (j in idx) {
  res <- simulate(post_b[j], post_g[j])
  if (!is.null(res)) curves[[length(curves) + 1]] <- res
}
curves_mat <- do.call(rbind, curves)
mean_pred <- colMeans(curves_mat)
lower_pred <- apply(curves_mat, 2, quantile, probs = 0.025)
upper_pred <- apply(curves_mat, 2, quantile, probs = 0.975)

# ── Plotting Styles ────────────────────────────────────────────────────────────
NAVY <- "#0D1B2A"; TEAL <- "#00A8CC"; ORANGE <- "#FF8C42"
MINT <- "#5CE1E6"; WHITE <- "#FFFFFF"; BGDARK <- "#132236"; BGLIGHT <- "#1C3F6E"

dark_theme <- theme_minimal() + theme(
  plot.background = element_rect(fill = NAVY, color = NA),
  panel.background = element_rect(fill = BGDARK, color = NA),
  panel.grid.major = element_line(color = BGLIGHT, linewidth = 0.5), 
  panel.grid.minor = element_line(color = BGLIGHT, linewidth = 0.25), 
  text = element_text(color = '#8FA8BE'),
  axis.text = element_text(color = '#8FA8BE', size = 9),
  axis.title = element_text(color = '#8FA8BE'),
  plot.title = element_text(color = WHITE, face = "bold", size = 13),
  plot.subtitle = element_text(color = WHITE),
  legend.background = element_rect(fill = BGDARK, color = NA),
  legend.text = element_text(color = WHITE),
  legend.title = element_text(color = WHITE)
)

# ── Fig 1: Trace Plots ─────────────────────────────────────────────────────────
df_trace_b <- data.frame(iter = 1:n_iter, value = chain_b, param = "β (Transmission Rate)")
df_trace_g <- data.frame(iter = 1:n_iter, value = chain_g, param = "γ (Recovery Rate)")

p_trace_b <- ggplot(df_trace_b, aes(x = iter, y = value)) +
  geom_line(color = TEAL, alpha = 0.88, linewidth = 0.45) + 
  geom_vline(xintercept = burn_in, color = WHITE, linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = mean(post_b), color = "yellow", linetype = "dotted", linewidth = 0.8) +
  labs(y = "β (Transmission Rate)", x = NULL, title = "MCMC Trace Plots — India COVID-19 Pre-Lockdown Phase") +
  dark_theme

p_trace_g <- ggplot(df_trace_g, aes(x = iter, y = value)) +
  geom_line(color = ORANGE, alpha = 0.88, linewidth = 0.45) + 
  geom_vline(xintercept = burn_in, color = WHITE, linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = mean(post_g), color = "yellow", linetype = "dotted", linewidth = 0.8) +
  labs(y = "γ (Recovery Rate)", x = "MCMC Iteration") +
  dark_theme

fig_trace <- p_trace_b / p_trace_g
print(fig_trace)

# ── Fig 2: Posteriors ──────────────────────────────────────────────────────────
plot_post <- function(data, title, col, xlab) {
  mn <- mean(data); lo <- quantile(data, 0.025); hi <- quantile(data, 0.975)
  df <- data.frame(val = data)
  
  ggplot(df, aes(x = val)) +
    geom_histogram(aes(y = after_stat(density)), bins = 55, fill = col, alpha = 0.83) +
    geom_vline(xintercept = mn, color = WHITE, linewidth = 1.2) + 
    geom_vline(xintercept = lo, color = "#aaaaaa", linewidth = 0.8, linetype = "dashed") +
    geom_vline(xintercept = hi, color = "#aaaaaa", linewidth = 0.8, linetype = "dashed") +
    labs(title = title, x = xlab, y = "Density") +
    dark_theme + theme(plot.title = element_text(color = col, size = 10.5))
}

p_post_b <- plot_post(post_b, 'β — Transmission Rate', TEAL, 'β')
p_post_g <- plot_post(post_g, 'γ — Recovery Rate', ORANGE, 'γ')
p_post_r <- plot_post(post_R0, 'R₀ = β / γ', MINT, 'R₀')

fig_post <- (p_post_b | p_post_g | p_post_r) + 
  plot_annotation(title = "Posterior Distributions — MCMC Estimates", theme = dark_theme)
print(fig_post)

# ── Fig 3: Observed vs Fitted ──────────────────────────────────────────────────
df_fit <- data.frame(
  day = 0:(T_len-1),
  obs = obs,
  mean_pred = mean_pred,
  lower_pred = lower_pred,
  upper_pred = upper_pred
)

date_labels <- c('Mar 14','Mar 17','Mar 20','Mar 23','Mar 26\n(Lockdown)','Mar 29','Apr 1')
tick_pos <- c(0, 3, 6, 9, 12, 15, 18)

fig_fit <- ggplot(df_fit, aes(x = day)) +
  geom_ribbon(aes(ymin = lower_pred, ymax = upper_pred), fill = TEAL, alpha = 0.28) +
  geom_line(aes(y = mean_pred), color = TEAL, linewidth = 1.2) + 
  geom_point(aes(y = obs), color = ORANGE, size = 3) + 
  geom_vline(xintercept = 11, color = WHITE, linetype = "dotted", linewidth = 0.8, alpha = 0.6) +
  annotate("text", x = 11.3, y = max(obs)*0.85, label = "Lockdown\nannounced\nMar 24", color = "#aaaaaa", size = 3, hjust = 0) +
  scale_x_continuous(breaks = tick_pos[tick_pos < T_len], labels = date_labels[tick_pos < T_len]) +
  labs(title = "India COVID-19 Pre-Lockdown: Observed Cases vs MCMC-Fitted SIR Model",
       subtitle = "(March 14–31, 2020 | 18 days | Real MoHFW data)",
       x = NULL, y = "Daily New Cases") +
  annotate("label", x = T_len * 0.60, y = max(obs) * 0.45, 
           label = sprintf("R₀ = %.2f\n[95%% CI: %.2f–%.2f]\n\nImplied infectious\nperiod: %.0f days\n\nData: MoHFW India / JHU CSSE", 
                           mean(post_R0), ci_R0[1], ci_R0[2], 1/mean(post_g)),
           color = MINT, fill = BGLIGHT, size = 3.5, label.padding = unit(0.5, "lines"), label.r = unit(0.2, "lines")) +
  dark_theme

print(fig_fit)

# ── Fig 4: Joint posterior (shows identifiability ridge) ───────────────────────
df_joint <- data.frame(
  beta = post_b[seq(1, length(post_b), by = 4)],
  gamma = post_g[seq(1, length(post_g), by = 4)],
  R0 = post_R0[seq(1, length(post_R0), by = 4)]
)

fig_joint <- ggplot(df_joint, aes(x = beta, y = gamma, color = R0)) +
  geom_point(alpha = 0.40, size = 1) +
  scale_color_viridis_c(option = "plasma", name = "R₀ = β/γ") +
  annotate("point", x = mean(post_b), y = mean(post_g), color = "white", shape = 8, size = 5) +
  labs(title = "Joint Posterior Distribution\nβ vs γ (colored by implied R₀)",
       subtitle = "— note the diagonal ridge: partial identifiability",
       x = "β (transmission rate)", y = "γ (recovery rate)") +
  dark_theme + theme(legend.position = "right")

print(fig_joint)

cat("\n✅ All plots generated and displayed!\n")