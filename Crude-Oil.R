# Title: How Does the U.S. Natural Gas Market React to Demand and Supply Shocks in the Crude Oil Market?
# Author: Behnam Fallah
# Date: 2025
# Notes:
#   - This is a cleaned & personalized rewrite of my original script,
#     keeping the overall structure (VAR → IRFs → cumulative effects → tables/plots)
#     I have improved naming, reproducibility, and readability.
#Used GPt-5 IN THE 0 & 1 part of the codes.
# ========================= 0) Setup & Libraries ===============================
rm(list = ls(all = TRUE))  # keep: your original workflow clears environment
set.seed(1234)

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(here)        # reproducible paths
  library(janitor)     # clean_names()
  library(vars)        # VAR estimation & IRFs
  library(Matrix)      # Cholesky, matrix ops
})

# Helper: concise error printing
fail <- function(...) stop(sprintf(...), call. = FALSE)

# ========================= 1) Data Ingest =====================================
# Personal touch: try multiple likely locations so future-you can run this anywhere.
DATA_PATHS <- c(
  here::here("data.txt"),
  here::here("data", "data.txt"),
  "data.txt"
)

existing <- DATA_PATHS[file.exists(DATA_PATHS)]
if (length(existing) == 0) fail("Could not find data.txt. Tried: %s", paste(DATA_PATHS, collapse = ", "))
DATA_FILE <- existing[[1]]

# Expecting a whitespace- or comma-delimited file with headers
raw_df <- tryCatch({
  readr::read_table2(DATA_FILE, col_types = cols(.default = col_double()))
}, error = function(e) {
  readr::read_csv(DATA_FILE, show_col_types = FALSE)
}) %>%
  janitor::clean_names()

# Minimal schema sanity checks (adjust to your real columns)
expected_cols <- c("prod", "ea", "po", "pg")
if (!all(expected_cols %in% names(raw_df))) {
  fail("Expected columns not found. Need: %s. Have: %s", 
       paste(expected_cols, collapse = ", "), paste(names(raw_df), collapse = ", "))
}

# Convert to ts object (monthly data assumed; adjust frequency/start as needed)
# Personalization: default start = c(1978, 1) to mirror your original timeline.
start_year <- 1978; start_period <- 1; freq <- 12
Y_ts <- ts(raw_df[, expected_cols], start = c(start_year, start_period), frequency = freq)
colnames(Y_ts) <- expected_cols

# ========================= 2) VAR Specification ===============================
# Choose lag length (keep your original p if known; otherwise select via AIC)
# Personal default: p = 12 for monthly dynamics, but also show data-driven choice.

p_default <- 12
lag_sel <- vars::VARselect(Y_ts, lag.max = 18, type = "const")
p_aic <- lag_sel$selection["AIC(n)"]
p <- ifelse(is.na(p_aic), p_default, as.integer(p_aic))

message(sprintf("Selected VAR lag length p = %d (AIC, fallback %d)", p, p_default))

# Estimate VAR
var_fit <- vars::VAR(Y_ts, p = p, type = "const")

# Core objects
Uhat <- residuals(var_fit)           # reduced-form residuals (T x k)
Sigma <- crossprod(Uhat) / nrow(Uhat) # residual covariance (k x k)
B <- chol(Sigma)                      # Cholesky factor (upper-tri); structural mapping
Ehat <- t(solve(B, t(Uhat)))         # structural shocks: U = B * e  ⇒  e = B^{-1} * U

# ========================= 3) Impulse Response Functions ======================
# Compute IRFs using vars::irf; keep horizon = 24 months (2 years) by default.
H <- 24
irf_all <- vars::irf(var_fit, impulse = colnames(Y_ts), response = colnames(Y_ts),
                     n.ahead = H, ortho = TRUE, boot = FALSE)

# Extract IRFs into tidy frames for plotting & later use
irf_list <- imap(irf_all$irf, ~ as_tibble(.x) %>% mutate(h = row_number() - 1, impulse = .y))
IRF_df <- bind_rows(irf_list) %>%
  pivot_longer(cols = all_of(expected_cols), names_to = "response", values_to = "irf") %>%
  relocate(impulse, response, h, irf)

# ========================= 4) Cumulative Effects on Oil Price (po) ============
# Your original computed cumulative contributions of prod, ea, po shocks to po.
# We reproduce that logic cleanly and transparently.

# Build IRF arrays for convenience (horizon H+1 because h=0..H)
get_irf_mat <- function(impulse_name) {
  IRF_df %>%
    filter(impulse == impulse_name) %>%
    arrange(h) %>%
    select(response, h, irf) %>%
    pivot_wider(names_from = response, values_from = irf) %>%
    as.matrix()
}

IRF_prod <- get_irf_mat("prod")   # rows: h, cols: prod/ea/po/pg
IRF_ea   <- get_irf_mat("ea")
IRF_po   <- get_irf_mat("po")
IRF_pg   <- get_irf_mat("pg")

# Cumulative effect of each shock on po; align with sample for which residuals exist
T_eff <- nrow(Uhat)   # = T - p

cum_effect <- function(IRF_mat, shock_series, target_col = "po") {
  # IRF_mat: (H+1) x k; shock_series: length T_eff of structural shocks for given impulse
  # We compute cumulative effect at time t as sum_{i=1..t} IRF[i, target]*e_t-i+1
  target_idx <- which(colnames(Y_ts) == target_col)
  yhat <- numeric(T_eff)
  for (i in seq_len(T_eff)) {
    # usable horizon limited by i
    h_use <- seq_len(min(i, nrow(IRF_mat)))
    yhat[i] <- sum(IRF_mat[h_use, target_idx] * rev(shock_series[seq_len(length(h_use))]))
  }
  yhat
}

yhat_prod <- cum_effect(IRF_prod, Ehat[, "prod"], target_col = "po")
yhat_ea   <- cum_effect(IRF_ea,   Ehat[, "ea"],   target_col = "po")
yhat_po   <- cum_effect(IRF_po,   Ehat[, "po"],   target_col = "po")

CumEffects_ts <- ts(cbind(yhat_prod, yhat_ea, yhat_po), start = start(Y_ts)[1:2] + c(0, p), frequency = freq)
colnames(CumEffects_ts) <- c("cum_prod_on_po", "cum_ea_on_po", "cum_po_on_po")

# ========================= 5) Plots (tidy) ====================================
plot_lines <- function(ts_obj, title, ylab = "", n_max = Inf) {
  df <- as_tibble(ts_obj) %>% mutate(t = zoo::as.yearmon(time(ts_obj)))
  if (is.finite(n_max)) df <- tail(df, n_max)
  df_long <- df %>% pivot_longer(-t, names_to = "series", values_to = "value")
  ggplot(df_long, aes(t, value, group = series)) +
    geom_line() +
    facet_wrap(~ series, scales = "free_y", ncol = 1) +
    labs(title = title, x = NULL, y = ylab) +
    theme_minimal(base_size = 12)
}

# Series levels
levels_plot <- plot_lines(Y_ts, title = "Series Levels: prod, ea, po, pg", ylab = "Level")
# IRFs to horizon H (po responses across impulses visible via filtering)
irf_plot_po <- IRF_df %>% filter(response == "po") %>%
  ggplot(aes(h, irf, color = impulse)) + geom_line() +
  labs(title = sprintf("Orthogonalized IRFs to horizon H = %d (response = po)", H), x = "Horizon (months)", y = "IRF") +
  theme_minimal(base_size = 12)

# Cumulative effects on po
cum_plot <- plot_lines(CumEffects_ts, title = "Cumulative Effect on Oil Price (po)", ylab = "Cumulative contribution")

# Print to viewer
print(levels_plot)
print(irf_plot_po)
print(cum_plot)

# ========================= 6) Tables & Summaries ==============================
# Example: contribution magnitudes at end-of-sample
end_contrib <- tail(as_tibble(CumEffects_ts), 1) %>%
  pivot_longer(everything(), names_to = "component", values_to = "value") %>%
  arrange(desc(abs(value)))

print(end_contrib)

# Optional: variance decomposition (FEVD)
fevd_obj <- vars::fevd(var_fit, n.ahead = H)
# Tidy FEVD for po
fevd_po <- fevd_obj$po %>%
  as_tibble() %>%
  mutate(h = row_number() - 1) %>%
  pivot_longer(cols = all_of(expected_cols), names_to = "impulse", values_to = "share") %>%
  arrange(h)

print(head(fevd_po, 10))

# ========================= 7) Save Outputs (optional) =========================
# ggsave(here::here("output", "levels_plot.png"), levels_plot, width = 7, height = 6, dpi = 300)
# ggsave(here::here("output", "irf_po_plot.png"), irf_plot_po, width = 7, height = 4, dpi = 300)
# ggsave(here::here("output", "cum_po_plot.png"), cum_plot, width = 7, height = 6, dpi = 300)

# ========================= 8) Takeaways (console) =============================
cat("\n--- EXECUTIVE SUMMARY (Behnam) ---\n")
cat(sprintf("VAR lag length (AIC): p = %d\n", p))
cat("IRFs computed to horizon:", H, "months\n")
cat("Displayed cumulative contributions of prod, ea, po shocks to oil price (po).\n")
cat("See 'end_contrib' table for magnitudes at sample end and 'fevd_po' for FEVD shares.\n")
