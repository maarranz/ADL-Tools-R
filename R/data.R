#' Simulated quarterly data with a known causal structure
#'
#' Three series, 1990 Q1 to 2019 Q4 (120 quarters), simulated from
#' \deqn{x_t = 0.5 x_{t-1} + e_{xt},}
#' \deqn{y_t = 1 + 0.5 y_{t-1} + 0.4 x_t + 0.3 x_{t-1} + e_{yt},}
#' \deqn{w_t = 0.6 w_{t-1} + e_{wt},}
#' with independent standard normal errors. So \eqn{x} Granger-causes
#' \eqn{y} with instantaneous causality and no feedback (type 5,
#' \eqn{x \Rightarrow y}); \eqn{w} is unrelated to both (type 1). The true
#' long-run multiplier of \eqn{x} on \eqn{y} is \eqn{(0.4 + 0.3)/(1 - 0.5) =
#' 1.4}.
#'
#' @format A quarterly `ts` matrix with 120 rows and columns `y`, `x`, `w`.
#' @source Simulated; the script is `data-raw/adl_sim.R` in the source
#'   repository.
#' @examples
#' adl(y ~ x, data = adl_sim, order = c(1, 1))
#' granger_test("y", "x", p = 2, data = adl_sim)
"adl_sim"
