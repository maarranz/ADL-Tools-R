# adltools 0.2.0

* Residual diagnostics: `bg_test()` (Breusch-Godfrey, LM and F forms),
  `arch_test()` (ARCH-LM), `reset_test()` (RESET, conventional or HAC) and
  `jb_test()` (Jarque-Bera), collected by `adl_diagnostics()`.
  Breusch-Godfrey and ARCH-LM always use OLS, as ADL3 explains.
* Stability: `chow_break_test()` (interaction form, HAC possible) and
  `chow_forecast_test()`, static (1-step, exact Omega with parameter
  uncertainty) or dynamic (N-step path, Omega from the psi-weights).
* `comfac_test()`: common factor restrictions of any order r (r * n
  restrictions), LR by concentrated least squares on the same sample; Wald
  form alpha_1 beta_0 + beta_1 = 0 for the ADL(1,1), as in ADL2.

# adltools 0.1.0

* `tidy()` and `glance()` methods for fits, tests, Granger results and
  dynamic effects (generics from the `generics` package).
* `adl_sim`: simulated quarterly data with a known causal structure.
* Vignette `granger`: Granger non-causality step by step.
* tsibble input: quarterly and monthly indexes are converted to zoo classes.
* `adl()` estimates ADL models by OLS from `order` + `first_lag`, an explicit
  `lags` list or a dynlm-style formula with `L()` terms; unlagged terms after
  `|`; data as ts, zoo/xts, data frame, tibble or tsibble.
* `adl_vcov()`: conventional or Newey-West HAC covariance, with the
  forecastingtools conventions.
* `adl_ic()`: AIC, BIC and Hannan-Quinn (EViews per-observation form) over
  a grid of lag orders on a common sample, with and without the
  contemporaneous regressors (`contemporaneous = "both"` by default).
* `summary()` of a fit shows z p-values by default (`dist = "t"` for the
  `summary.lm()` values), takes `cov_type` for HAC standard errors, and ends
  with a Wald test that all slopes are zero. Coefficient tables label each
  parameter in the chapters' notation (`alpha_1`, `beta_0`, `beta_2,1`, ...).
* `adl_dynamics()`: dynamic and cumulative multipliers (also as a percentage
  of the long run), mean and median lag, and the inverse roots of C(L), with
  `print()` and `plot()` methods. Same definitions as adl-theory.
* Inference uses the asymptotic normal and chi-square distributions by
  default, with either covariance; `dist = "t"` or `"F"` gives the
  EViews/`summary.lm()` versions.
* `granger_test()` and `granger_test_general()`: Granger non-causality
  tests by Sargent's direct procedure (four models on a common sample,
  chi-square and F statistics) and the eight-way classification of causal
  relations.
* Inference: `adl_coefficients()` (coefficient table), `adl_joint_test()`
  (Wald tests on coefficients, whole variables or `R b = r`) and
  `adl_long_run()` (long-run multipliers with delta-method standard errors).
* Package skeleton and design document (`DESIGN.md`).
* Design notes updated after the ADL1–ADL4 revisions: several regressors and
  per-variable lag sets, contemporaneous-term rule, generalized COMFAC test,
  exogeneity tests, static/dynamic Chow forecast test, Gets, revised roadmap.
