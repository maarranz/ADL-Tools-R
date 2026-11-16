# adltools (development version)

* `adl()` estimates ADL models by OLS from `order` + `first_lag`, an explicit
  `lags` list or a dynlm-style formula with `L()` terms; unlagged terms after
  `|`; data as ts, zoo/xts, data frame, tibble or tsibble.
* `adl_vcov()`: conventional or Newey-West HAC covariance, with the
  forecastingtools conventions.
* Package skeleton and design document (`DESIGN.md`).
* Design notes updated after the ADL1–ADL4 revisions: several regressors and
  per-variable lag sets, contemporaneous-term rule, generalized COMFAC test,
  exogeneity tests, static/dynamic Chow forecast test, Gets, revised roadmap.
