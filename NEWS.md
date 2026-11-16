# adltools (development version)

* `adl()` estimates ADL models by OLS from `order` + `first_lag`, an explicit
  `lags` list or a dynlm-style formula with `L()` terms; unlagged terms after
  `|`; data as ts, zoo/xts, data frame, tibble or tsibble.
* `adl_vcov()`: conventional or Newey-West HAC covariance, with the
  forecastingtools conventions.
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
