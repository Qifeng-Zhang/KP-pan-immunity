# ============================================================
# 04_Chemostat_parameter_scan.R
# Model   : dN/dt = D * (N_in - N) - (1/Y) * (rate_S * S + rate_R * R)
#           dS/dt = rate_S * S - phi * S * P - D * S
#           dR/dt = rate_R * R - D * R
#           dP/dt = beta * phi * S * P - lambda * R * P - D * P
#           cost = 1 - mu_max_R / mu_max_S ;  mu_max_R = mu_max_S * (1 - cost)
# Tool    : deSolve, foreach + doParallel
# Input   : parameters
# Output  : chemostat_scan_results.csv
#             cost, lambda, prop_S
# Manuscript: Main text Fig. 4 / Supplementary Note 5
# ============================================================

library(deSolve)
library(tidyverse)
library(foreach)
library(doParallel)

OUT_DIR <- "/path/to/model_output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Chemostat model
# -----------------------------
monod_chemostat <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    N <- max(N, 0); S <- max(S, 0); R <- max(R, 0); P <- max(P, 0)

    rate_S <- mu_max_S * (N / (Ks + N))
    rate_R <- mu_max_R * (N / (Ks + N))
    lysis_flux <- phi * S * P

    dN <- D * (N_in - N) - (1 / Y) * (rate_S * S + rate_R * R)
    dS <- rate_S * S - lysis_flux - D * S
    dR <- rate_R * R - D * R
    dP <- beta * lysis_flux - lambda * R * P - D * P

    list(c(dN, dS, dR, dP))
  })
}

# -----------------------------
# 2. Fixed parameters and scan grid
# -----------------------------
fixed_mu_S <- 0.57

parms_base <- c(
  mu_max_S = fixed_mu_S,
  Ks    = 0.29,
  Y     = 1.0,
  phi   = 0.40,
  beta  = 100,
  D     = 0.1,
  N_in  = 1.0,
  mu_max_R = NA,
  lambda   = NA
)

RESOLUTION <- 100
COST_RANGE   <- c(-0.1, 0.6)
LAMBDA_RANGE <- c(0, 40)

T_END     <- 240
T_MEASURE <- 120

grid <- expand.grid(
  cost   = seq(COST_RANGE[1],   COST_RANGE[2],   length.out = RESOLUTION),
  lambda = seq(LAMBDA_RANGE[1], LAMBDA_RANGE[2], length.out = RESOLUTION)
)

# -----------------------------
# 3. Parallel scan
# -----------------------------
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(
  i = 1:nrow(grid),
  .combine = rbind,
  .packages = c("deSolve", "dplyr"),
  .export = c("monod_chemostat", "parms_base", "fixed_mu_S")
) %dopar% {

  curr_parms <- parms_base
  curr_parms["mu_max_R"] <- fixed_mu_S * (1 - grid$cost[i])
  curr_parms["lambda"]   <- grid$lambda[i]

  out <- tryCatch({
    ode(y = c(N = 0.5, S = 0.030, R = 0.030, P = 1),
        times = seq(0, T_END, by = 2),
        func = monod_chemostat, parms = curr_parms,
        method = "lsoda", maxsteps = 200000, atol = 1e-7, rtol = 1e-7)
  }, error = function(e) return(NULL))

  prop_S <- NA
  if (!is.null(out)) {
    steady <- subset(as.data.frame(out), time >= T_MEASURE)
    mean_S <- mean(steady$S, na.rm = TRUE)
    mean_R <- mean(steady$R, na.rm = TRUE)
    prop_S <- mean_S / (mean_S + mean_R)
  }

  data.frame(cost = grid$cost[i], lambda = grid$lambda[i], prop_S = prop_S)
}

stopCluster(cl)

# -----------------------------
# 4. Export results
# -----------------------------
write.csv(results, file.path(OUT_DIR, "chemostat_scan_results.csv"), row.names = FALSE)
