# =============================================================================
# cosinor_analysis.R
#
# Utilities for fitting a cosinor (cosine) model to unevenly spaced time
# series data with an unknown (free) period.
#
# Model:  y = A * cos(2π * x / T + φ) + C + ε
#
#   T   – period
#   A   – amplitude
#   φ   – phase shift (radians)
#   C   – MESOR (midline estimating statistic of rhythm; the mean level)
#   ε   – zero-mean Gaussian noise
#
# Functions
#   generate_timeseries()        Simulate an irregularly sampled cosine signal
#   fit_cosinor_free_period()    Estimate T, A, φ, C from data
#   print_cosinor_summary()      Pretty-print the fit results
# =============================================================================


# -----------------------------------------------------------------------------
# generate_timeseries
#
# Simulate an unevenly sampled cosine time series with additive Gaussian noise.
#
# Arguments:
#   n         Number of observations to attempt before trimming (default 500)
#   T_true    True period of the cosine signal             (default 24)
#   A_true    True amplitude                               (default 10)
#   phi_true  True phase shift in radians                  (default 0.5)
#   C_true    True MESOR (baseline / mean level)           (default 60)
#   noise_sd  Standard deviation of the additive noise     (default 3)
#   Tmax      Maximum observation time; observations beyond
#             this value are discarded                     (default 72)
#   seed      Random seed for reproducibility              (default 1)
#
# Returns:
#   A data.frame with columns:
#     x       – observation times (irregularly spaced, in [0, Tmax])
#     y       – noisy observations
#     y_true  – noise-free signal values
# -----------------------------------------------------------------------------
generate_timeseries <- function(n       = 500,
                                T_true  = 24,
                                A_true  = 10,
                                phi_true = 0.5,
                                C_true  = 60,
                                noise_sd = 3,
                                Tmax    = 72,
                                seed    = 1) {
  set.seed(seed)

  # Draw inter-arrival times from Exp(n / Tmax) so that the expected total
  # span is approximately Tmax.
  dt <- rexp(n, rate = n / Tmax)
  x  <- cumsum(dt)

  # Keep only observations within [0, Tmax].
  x <- x[x <= Tmax]
  n <- length(x)

  # Compute the true (noise-free) signal and add Gaussian noise.
  y_true <- A_true * cos(2 * pi * x / T_true + phi_true) + C_true
  y      <- y_true + rnorm(n, sd = noise_sd)

  data.frame(x = x, y = y, y_true = y_true)
}


# -----------------------------------------------------------------------------
# fit_cosinor_free_period
#
# Fit a cosinor model with an unknown period using a two-stage optimisation:
#   1. Coarse grid search over T_range to locate the approximate minimum RSS.
#   2. Fine-scale optimisation (golden-section search) in a narrow interval
#      around the best grid point.
#
# For a fixed candidate period T0, the model is linear in its remaining
# parameters (C, β_cos, β_sin), so ordinary least squares is used via
# lm.fit() at each function evaluation.
#
# Arguments:
#   x                Numeric vector of observation times
#   y                Numeric vector of observed values (same length as x)
#   T_range          Two-element numeric vector [T_min, T_max] defining the
#                    search range for the period             (default c(20, 28))
#   T_grid_step      Step size for the coarse grid search   (default 10/60)
#   refine_half_width Half-width of the refinement interval around the coarse
#                    optimum                                 (default 0.5)
#   tol              Convergence tolerance for optimize()   (default 1e-5)
#
# Returns:
#   A named list with the following elements:
#     T            – estimated period
#     A            – estimated amplitude  (>= 0)
#     phi          – estimated phase shift (radians, in (-π, π])
#     C            – estimated MESOR
#     coefficients – named vector (intercept, beta_cos, beta_sin)
#     fitted.values – fitted values at each observation
#     residuals    – residuals (y - fitted)
#     RSS          – residual sum of squares
#     R_squared    – coefficient of determination
#     T_grid       – period grid used in the coarse search
#     RSS_grid     – RSS evaluated at each grid point
#     method       – character string describing the estimation approach
# -----------------------------------------------------------------------------
fit_cosinor_free_period <- function(x,
                                    y,
                                    T_range          = c(20, 28),
                                    T_grid_step      = 10/60,
                                    refine_half_width = 0.5,
                                    tol              = 1e-5) {

  # ---- Input validation ----
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors.")
  }
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }

  # Remove non-finite observations.
  sel <- is.finite(x) & is.finite(y)
  x   <- x[sel]
  y   <- y[sel]

  if (length(x) < 4) {
    stop("At least 4 finite data points are required.")
  }

  # ---- RSS as a function of the period T0 ----
  # For fixed T0, regress y on the design matrix [1, cos(2πx/T0), sin(2πx/T0)]
  # using OLS, and return the residual sum of squares.
  rss_fun <- function(T0) {
    X   <- cbind(1,
                 cos(2 * pi * x / T0),
                 sin(2 * pi * x / T0))
    fit <- lm.fit(x = X, y = y)
    sum(fit$residuals^2)
  }

  # ---- Stage 1: Coarse grid search ----
  T_grid   <- seq(T_range[1], T_range[2], by = T_grid_step)
  rss_grid <- vapply(T_grid, rss_fun, numeric(1))
  T0       <- T_grid[which.min(rss_grid)]   # best grid candidate

  # ---- Stage 2: Fine-scale optimisation ----
  lower <- max(T_range[1], T0 - refine_half_width)
  upper <- min(T_range[2], T0 + refine_half_width)

  opt   <- optimize(rss_fun, interval = c(lower, upper), tol = tol)
  T_opt <- opt$minimum

  # ---- Final OLS fit at the optimal period ----
  X   <- cbind(1,
               cos(2 * pi * x / T_opt),
               sin(2 * pi * x / T_opt))
  fit <- lm.fit(x = X, y = y)

  C     <- fit$coefficients[1]   # intercept  = MESOR
  beta1 <- fit$coefficients[2]   # coefficient of cos term
  beta2 <- fit$coefficients[3]   # coefficient of sin term

  # Convert linear coefficients to amplitude and phase:
  #   A * cos(θ + φ) = A*cos(φ)*cos(θ) - A*sin(φ)*sin(θ)
  #   so β1 = A*cos(φ),  β2 = -A*sin(φ)
  A   <- sqrt(beta1^2 + beta2^2)
  phi <- atan2(-beta2, beta1)

  # Goodness-of-fit statistics
  y_hat     <- as.vector(X %*% fit$coefficients)
  residuals <- y - y_hat
  rss       <- sum(residuals^2)
  tss       <- sum((y - mean(y))^2)
  r_squared <- 1 - rss / tss

  list(
    T             = T_opt,
    A             = A,
    phi           = phi,
    C             = C,
    coefficients  = c(intercept = C, beta_cos = beta1, beta_sin = beta2),
    fitted.values = y_hat,
    residuals     = residuals,
    RSS           = rss,
    R_squared     = r_squared,
    T_grid        = T_grid,
    RSS_grid      = rss_grid,
    method        = "coarse grid + optimize + linear regression"
  )
}


# -----------------------------------------------------------------------------
# print_cosinor_summary
#
# Display a formatted summary of a cosinor fit produced by
# fit_cosinor_free_period().
#
# Arguments:
#   fit     List returned by fit_cosinor_free_period()
#   digits  Number of decimal places for printed values   (default 3)
#
# Returns:
#   NULL (invisibly). Called for its side-effect of printing to the console.
# -----------------------------------------------------------------------------
print_cosinor_summary <- function(fit, digits = 3) {

  if (is.null(fit$T)) {
    stop("Invalid fit object: element 'T' is missing.")
  }

  # Extract estimates.
  T   <- fit$T
  A   <- fit$A
  phi <- fit$phi
  C   <- fit$C
  R2  <- fit$R_squared

  # Convert phase to peak time:
  #   Maximum of cos(2πx/T + φ) occurs when 2πx/T + φ = 0,
  #   i.e., x_peak = -φ / (2π) * T.
  # Wrap into [0, T) so the result is always interpretable as a time
  # within the first cycle.
  t_peak <- ((-phi / (2 * pi)) * T) %% T

  # ---- Print ----
  cat("=====================================\n")
  cat(" Cosinor Fit Summary (Free Period)\n")
  cat("=====================================\n\n")

  cat("[Model]\n")
  cat(" y = A * cos(2*pi*x / T + phi) + C\n\n")

  cat("[Estimated Parameters]\n")
  cat(sprintf(" Period    (T)   : %.*f\n",        digits, T))
  cat(sprintf(" Amplitude (A)   : %.*f\n",        digits, A))
  cat(sprintf(" Phase     (phi) : %.*f  rad\n",   digits, phi))
  cat(sprintf(" MESOR     (C)   : %.*f\n",        digits, C))
  cat(sprintf(" Peak time       : %.*f  (same units as x)\n\n", digits, t_peak))

  cat("[Goodness of Fit]\n")
  cat(sprintf(" R-squared : %.*f\n",  digits, R2))
  cat(sprintf(" RSS       : %.*f\n\n", digits, fit$RSS))

  cat("[Interpretation]\n")
  cat(sprintf(" T   ~ %.*f : estimated period\n",                   digits, T))
  cat(sprintf(" A   ~ %.*f : oscillation magnitude\n",              digits, A))
  cat(sprintf(" phi ~ %.*f : phase shift (rad)\n",                  digits, phi))
  cat(sprintf(" t_peak ~ %.*f : time at which the peak occurs\n",   digits, t_peak))
  cat(sprintf(" C   ~ %.*f : mean level (MESOR)\n",                 digits, C))
  cat(sprintf(" R^2 ~ %.*f : proportion of variance explained\n",   digits, R2))

  cat("\n=====================================\n")

  invisible(NULL)
}
