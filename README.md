# adltools (R)

[![R-CMD-check](https://github.com/maarranz/ADL-Tools-R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/maarranz/ADL-Tools-R/actions/workflows/R-CMD-check.yaml)

Teaching-oriented R tools for **autoregressive distributed lag (ADL)** models, developed for **Applied Time Series Econometrics (ATSE)** at **Universidad Carlos III de Madrid**.

> **Status: version 0.5.0.** Formulation, estimation, inference, the ECM form, dynamic effects with standard errors, the J-test for non-nested models, Granger non-causality (with automatic lag selection), residual diagnostics, CUSUM, CUSUM-of-squares, Hansen and Chow tests, the generalized COMFAC test and exogeneity tests (weak, strong, super), outlier and break detection, and general-to-specific lag selection. The design decisions are recorded in [`DESIGN.md`](DESIGN.md).

## Scope

`adltools` brings together, under one consistent interface, the tools needed to formulate, estimate and test ADL models as they are taught in the ATSE chapters on dynamic regression (ADL I–IV). Several of these tests are scattered across R packages, and some are not available elsewhere in the form used in the course.

| Task | Functions |
|:--|:--|
| Formulation and estimation | `adl()`: lag orders with `order` and `first_lag`, explicit lag sets, or `dynlm`-style formulas; several regressors; unlagged dummies after `\|`; `ts`, `zoo`/`xts`, data frames, tibbles and tsibbles |
| Inference | `summary()`, `adl_coefficients()` (parameters in the chapters' α/β notation), `adl_joint_test()`, `adl_long_run()`; conventional or Newey–West (HAC) covariance |
| Error correction form | `adl_ecm()`: the fitted ADL rewritten as an ECM (exact standard errors), with the speed of adjustment, the long-run multipliers and the equilibrium error |
| Dynamic effects | `adl_dynamics()`: dynamic and cumulative multipliers, mean and median lag, with standard errors; inverse roots |
| Lag orders | `adl_ic()`: AIC, BIC and Hannan–Quinn on a common sample, with and without the contemporaneous regressor |
| General-to-specific | `adl_gets()`: from ADL(p, p) down, joint Wald tests of the highest lags, diagnostics at every step |
| Granger non-causality | `granger_test()`, `granger_test_general()`: Sargent's direct procedure and the eight causal types; orders by BIC/AIC/HQ or Gets with `select` |
| Diagnostics | `adl_diagnostics()`, `bg_test()` (Breusch–Godfrey), `arch_test()` (ARCH-LM), `reset_test()` (RESET), `jb_test()` (Jarque–Bera) |
| Stability | `cusum_test()`, `cusumsq_test()` (with plots), `hansen_test()` (per parameter, no break date), `chow_break_test()`, `chow_forecast_test()` (static 1-step or dynamic N-step) |
| Non-nested models | `j_test()`: the J-test in both directions, with the four-outcome reading of the chapters; `adl` or `dynlm` fits |
| Common factors | `comfac_test()`: COMFAC of any order (LR), and the ADL(1,1) Wald form of the chapters |
| Outliers and breaks | `adl_outliers()` (impulse and step dummies; own search or `gets::isat()`) |
| Exogeneity | `wu_hausman_test()` (predeterminedness), `weak_exogeneity_test()` (weak, and strong with `strong = TRUE`), `super_exogeneity_test()` (invariance to breaks in the marginal process) |
| Existing fits | `as_adl()` converts a `dynlm` fit |
| Output | `print()`, `plot()`, `tidy()`, `glance()` |

**Scope.** The tests assume stationary variables, or, where the chapters say so, cointegrated ones. The package does not test for unit roots or cointegration and provides no bounds test or cointegrating-regression ECM (`adl_ecm()` only rewrites a fitted ADL): unit roots belong to `unitroottools`, and cointegration to a planned companion package. With integrated series, a finite long-run multiplier or a conventional Wald p-value does not establish a long-run relation.

Inference uses the asymptotic normal and chi-square distributions by default: the regressions are dynamic, so exact t and F results do not apply. `dist = "t"` or `dist = "F"` reproduces the figures of EViews and `summary.lm()`.

## Example

```r
library(adltools)

# Three ways to write the same ADL(1, 1)
fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
fit <- adl(y ~ x, data = adl_sim, lags = list(y = 1, x = 0:1))
fit <- adl(y ~ L(y, 1) + L(x, 0:1), data = adl_sim)

summary(fit)                         # alpha/beta labels, z p-values
summary(fit, cov_type = "HAC")       # Newey-West standard errors
adl_long_run(fit)                    # long-run multiplier, delta method
adl_ecm(fit)                         # the same model in ECM form
adl_dynamics(fit)                    # multipliers, mean and median lag, with SEs
adl_joint_test(fit, variables = "x") # Wald test: all lags of x are zero

# The dead-start model: x enters from lag 1
adl(y ~ x, data = adl_sim, order = c(1, 1), first_lag = 1)

# Choosing lag orders on a common sample
adl_ic(y ~ x, data = adl_sim, max_order = c(4, 4))

# Granger non-causality, both directions and instantaneous
granger_test("y", "x", p = 2, data = adl_sim)

# Diagnostics, stability and common factors
adl_diagnostics(fit, bg_order = 4)
chow_break_test(fit, c(2005, 1))
chow_forecast_test(fit, c(2016, 1), type = "dynamic")
comfac_test(fit)
```

Two vignettes: a complete analysis following the ADL chapters, `vignette("workflow", package = "adltools")`, and the Granger procedure, `vignette("granger", package = "adltools")`.

## Related packages

`adltools` is part of a suite of teaching packages for ATSE:

| Package | Repository | Covers |
|:--|:--|:--|
| arimatools | [ARIMA-Tools-R](https://github.com/maarranz/ARIMA-Tools-R) | ARIMA estimation, diagnostics, outliers and dynamic forecasting |
| forecastingtools | [Forecasting-Tools-R](https://github.com/maarranz/Forecasting-Tools-R) | Static forecasting, forecast accuracy, evaluation and comparison (including ADL forecasts) |
| unitroottools | [UnitRoot-Tools-R](https://github.com/maarranz/UnitRoot-Tools-R) | Unit-root and stationarity testing |
| adltools | ADL-Tools-R | ADL models (this package) |

## Installation

```r
# install.packages("remotes")
remotes::install_github("maarranz/ADL-Tools-R", build_vignettes = TRUE)
library(adltools)
```

Requires R ≥ 4.2.

## Author

**Miguel A. Arranz**
Universidad Carlos III de Madrid
<maarranz@eco.uc3m.es>
ORCID: [0000-0002-5951-2284](https://orcid.org/0000-0002-5951-2284)

## Licence

MIT © Miguel A. Arranz
