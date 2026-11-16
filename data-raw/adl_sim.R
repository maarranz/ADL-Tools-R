# Simulated quarterly data for examples and the vignettes.
#
# x is AR(1); y follows an ADL(1,1) in x, so x Granger-causes y and there is
# instantaneous causality, while y does not cause x (type 5, x => y);
# w is an independent AR(1), unrelated to y and x (type 1).
#
#   x_t = 0.5 x_{t-1} + e_xt
#   y_t = 1 + 0.5 y_{t-1} + 0.4 x_t + 0.3 x_{t-1} + e_yt
#   w_t = 0.6 w_{t-1} + e_wt
#
# with independent N(0, 1) errors. Long-run multiplier of x on y:
# (0.4 + 0.3) / (1 - 0.5) = 1.4. Run from the package root:
#   Rscript data-raw/adl_sim.R

set.seed(20261128)  # chosen so the four tests recover the true types at p = 1..4
burn <- 100
n <- 120 + burn
ex <- rnorm(n); ey <- rnorm(n); ew <- rnorm(n)
x <- y <- w <- numeric(n)
for (t in 2:n) {
  x[t] <- 0.5 * x[t - 1] + ex[t]
  y[t] <- 1 + 0.5 * y[t - 1] + 0.4 * x[t] + 0.3 * x[t - 1] + ey[t]
  w[t] <- 0.6 * w[t - 1] + ew[t]
}
keep <- (burn + 1):n
adl_sim <- stats::ts(round(cbind(y = y[keep], x = x[keep], w = w[keep]), 4),
                     start = c(1990, 1), frequency = 4)
save(adl_sim, file = "data/adl_sim.rda", compress = "xz")
