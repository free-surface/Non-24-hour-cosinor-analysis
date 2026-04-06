# =============================================================================
# cosinor_example.R
#
# Demonstration script for the cosinor analysis workflow:
#   1. Simulate an irregularly sampled cosine time series.
#   2. Fit a cosinor model with a free (unknown) period.
#   3. Plot the raw data together with the fitted curve.
#   4. Print a formatted summary of the estimated parameters.
#
# Dependencies: cosinor_analysis.R  (must be sourced before running this file)
# =============================================================================

# source("cosinor_analysis.R")   # uncomment if functions are not yet loaded


# -----------------------------------------------------------------------------
# 1. Simulate data
#
# Generate ~288 irregularly spaced observations over 72 hours.
# True model: y = 8 * cos(2π * x / 25.3 + 0.8) + 70 + N(0, 3.5²)
# -----------------------------------------------------------------------------
dat <- generate_timeseries(
  n        = 288,    # target number of observations before trimming
  T_true   = 25.3,  # true period (hours)
  A_true   = 8,     # true amplitude
  phi_true = 0.8,   # true phase shift (radians)
  C_true   = 70,    # true MESOR (mean level)
  noise_sd = 3.5,   # standard deviation of additive Gaussian noise
  Tmax     = 72,    # observation window (hours); points beyond this are dropped
  seed     = 12     # random seed for reproducibility
)


# -----------------------------------------------------------------------------
# 2. Fit the cosinor model
#
# Period is treated as unknown and searched over [20, 28] hours.
# A coarse grid search identifies the approximate optimum, which is then
# refined with a golden-section search (optimize()) to tolerance 1e-5.
# -----------------------------------------------------------------------------
fit <- fit_cosinor_free_period(
  x       = dat$x,
  y       = dat$y,
  T_range = c(20, 28),   # plausible period search range (hours)
  tol     = 1e-5         # convergence tolerance for the fine-scale optimiser
)


# -----------------------------------------------------------------------------
# 3. Plot: raw observations + fitted curve
# -----------------------------------------------------------------------------

# -- Raw observations --
plot(dat$x, dat$y,
     pch  = 16,            # filled circles
     col  = "gray60",      # light gray to keep fitted curve visually prominent
     xlab = "Time (hours)",
     ylab = "y",
     main = "Cosinor Fit (Free Period)")

# -- Fitted cosine curve --
# Evaluated as a continuous function over the observed time range so the curve
# looks smooth regardless of the (irregular) sampling grid.
curve(
  fit$A * cos(2 * pi * x / fit$T + fit$phi) + fit$C,
  from = min(dat$x),
  to   = max(dat$x),
  add  = TRUE,       # overlay on the existing plot
  col  = "firebrick",
  lwd  = 2,
  n    = 1000        # number of points used to draw the curve; increase for
                     # finer resolution if the period is very short
)

# -- Optional: add the true (noise-free) signal for visual comparison --
# curve(
#   8 * cos(2 * pi * x / 25.3 + 0.8) + 70,
#   from = min(dat$x), to = max(dat$x),
#   add = TRUE, col = "steelblue", lwd = 1.5, lty = 2, n = 1000
# )

# -- Legend --
legend("topright",
       legend = c("Observed", "Fitted"),
       pch    = c(16, NA),
       lty    = c(NA, 1),
       lwd    = c(NA, 2),
       col    = c("gray60", "firebrick"),
       bty    = "n")


# -----------------------------------------------------------------------------
# 4. Print parameter estimates and goodness-of-fit statistics
# -----------------------------------------------------------------------------
print_cosinor_summary(fit)