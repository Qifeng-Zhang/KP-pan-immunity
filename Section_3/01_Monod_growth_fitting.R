# ============================================================
# 01_Monod_growth_fitting.R
# Model   : dB/dt = mu_max * N / (Ks + N) * B ;  dN/dt = -dB/dt
#           with N0 = ODmax - OD0 and B0 = OD0
# Tool    : deSolve, FME
# Input   : growth curves
# Output  : monod_parameters.csv 
#           monod_group_summary.csv
#           monod_fitted_curves.csv
#           monod_observed_points.csv
# ============================================================

library(deSolve)
library(FME)
library(tidyverse)

OUT_DIR <- "/path/to/model_output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Growth-curve data
# -----------------------------
SAMPLE_INTERVAL <- 1.5
times_exp <- seq(0, 12, by = SAMPLE_INTERVAL)

OD_data <- data.frame(
  Time = times_exp,
  S1 = c(0.060, 0.182, 0.287, 0.406, 0.623, 1.011, 1.099, 1.093, 1.137),
  S2 = c(0.046, 0.125, 0.132, 0.382, 0.554, 0.926, 1.081, 1.076, 1.061),
  S3 = c(0.052, 0.156, 0.202, 0.507, 0.666, 1.063, 1.088, 1.086, 1.083),
  S4 = c(0.056, 0.096, 0.126, 0.242, 0.452, 0.758, 0.922, 1.042, 1.049),
  R1 = c(0.069, 0.112, 0.129, 0.211, 0.416, 0.551, 0.733, 0.951, 1.017),
  R2 = c(0.065, 0.065, 0.111, 0.232, 0.405, 0.501, 0.806, 0.901, 1.008),
  R3 = c(0.073, 0.096, 0.183, 0.276, 0.382, 0.503, 0.638, 0.877, 1.030)
)

long_data <- OD_data %>%
  pivot_longer(cols = -Time, names_to = "Strain", values_to = "OD")

# -----------------------------
# 2. Monod model and fitting
# -----------------------------
monod_model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    N <- max(N, 0); B <- max(B, 0)
    growth_rate <- mu_max * (N / (Ks + N)) * B
    list(c(-growth_rate, growth_rate))
  })
}

fit_strain <- function(strain_df) {
  B_init <- strain_df$OD[1]
  B_max  <- max(strain_df$OD)
  state_init <- c(N = B_max - B_init, B = B_init)

  cost_func <- function(pars) {
    out <- ode(y = state_init, times = strain_df$Time,
               func = monod_model, parms = pars)
    return(out[, "B"] - strain_df$OD)
  }

  fit <- modFit(f = cost_func, p = c(mu_max = 0.5, Ks = 0.2), lower = c(0, 0))

  times_smooth <- seq(0, 12, length.out = 100)
  out_smooth <- ode(y = state_init, times = times_smooth,
                    func = monod_model, parms = fit$par)
  pred_df <- as.data.frame(out_smooth) %>%
    select(time, B) %>%
    rename(Time = time, OD_pred = B)

  return(list(params = fit$par, pred = pred_df))
}

results_params <- data.frame()
results_pred   <- data.frame()

for (s in unique(long_data$Strain)) {
  fit_res <- fit_strain(long_data %>% filter(Strain == s))
  results_params <- rbind(
    results_params,
    data.frame(Strain = s,
               mu_max = fit_res$params["mu_max"],
               Ks     = fit_res$params["Ks"])
  )
  results_pred <- rbind(results_pred, transform(fit_res$pred, Strain = s))
}

summary_stats <- results_params %>%
  mutate(Group = ifelse(grepl("S", Strain), "Sensitive", "Resistant")) %>%
  group_by(Group) %>%
  summarise(Mean_mu = mean(mu_max), SD_mu = sd(mu_max), Mean_Ks = mean(Ks))

print(summary_stats)

# -----------------------------
# 3. Export results
# -----------------------------
write.csv(results_params, file.path(OUT_DIR, "monod_parameters.csv"), row.names = FALSE)
write.csv(summary_stats,  file.path(OUT_DIR, "monod_group_summary.csv"), row.names = FALSE)

write.csv(
  long_data %>% pivot_wider(names_from = Strain, values_from = OD) %>% arrange(Time),
  file.path(OUT_DIR, "monod_observed_points.csv"), row.names = FALSE
)
write.csv(
  results_pred %>% pivot_wider(names_from = Strain, values_from = OD_pred) %>% arrange(Time),
  file.path(OUT_DIR, "monod_fitted_curves.csv"), row.names = FALSE
)
