# ============================================================
# 03_Model_2_12h_simulation.R
# Model   : dN/dt = -(1/Y) * (rate_S * S + rate_R * R)
#           dS/dt = rate_S * S - phi * S * P
#           dR/dt = rate_R * R
#           dP/dt = beta * phi * S * P - lambda * R * P
#           rate_i = mu_max_i * N / (Ks_i + N)
# Tool    : deSolve
# Input   : parameters
# Output  : model2_12h_timeseries.csv
# ============================================================

library(deSolve)
library(tidyverse)

OUT_DIR <- "/path/to/model_output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Parameters
# -----------------------------
parms <- c(
  mu_max_S = 0.57,
  mu_max_R = 0.35,
  Ks_S     = 0.29,
  Ks_R     = 0.29,
  Y        = 1.0,
  phi      = 0.40,
  beta     = 100,
  lambda   = 20
)

# -----------------------------
# 2. Model definition
# -----------------------------
model_2 <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    N <- max(N, 0); S <- max(S, 0); R <- max(R, 0); P <- max(P, 0)

    rate_S <- mu_max_S * (N / (Ks_S + N))
    rate_R <- mu_max_R * (N / (Ks_R + N))

    dN <- -(1 / Y) * (rate_S * S + rate_R * R)
    dS <- rate_S * S - phi * S * P
    dR <- rate_R * R
    dP <- beta * phi * S * P - lambda * R * P

    list(c(dN, dS, dR, dP))
  })
}

# -----------------------------
# 3. Simulation
# -----------------------------
yini  <- c(N = 1, S = 0.03, R = 0.03, P = 1.00)
times <- seq(0, 12, by = 0.1)

out <- ode(y = yini, times = times, func = model_2, parms = parms)

# -----------------------------
# 4. Derived quantities and export
# -----------------------------
df <- as.data.frame(out) %>%
  mutate(
    Total_Bio = S + R,
    Rel_S = S / Total_Bio,
    Rel_R = R / Total_Bio
  )

write.csv(df, file.path(OUT_DIR, "model2_12h_timeseries.csv"), row.names = FALSE)
