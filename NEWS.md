# adltools 0.2.1.9000

* New vignette `workflow`: one ADL analysis from lag selection to the long
  run and the dynamics, following the ADL chapters.

## Fixes from the second audit

* `weak_exogeneity_test()`: a control or regressor whose name clashed with a
  column of the marginal regression (e.g. a control called `ec`) silently
  replaced the equilibrium error, so the reported test was on the wrong
  coefficient. The test's own columns now keep their names and others are
  renamed.
* Dates must be consecutive months, quarters or years *and* sit at the same
  point of each period (same day, or month end). Before, 1 Jan, 25 Feb and
  1 Mar were accepted as three consecutive months. The gap message now
  points at the actual gap.
* Gaps are checked only within the observations used (the estimation window
  and its initial lags), so a gap elsewhere in the data no longer prevents
  estimating on a regular stretch with `start`/`end`.
* `adl_dynamics()` with a `max_horizon` shorter than the distributed lag
  failed with an internal error; it now warns and stops at `max_horizon`.

# adltools 0.2.1

Fixes from the external audit of 2026-11-16 (medium priority):

* `adl_ic()` drops collinear or saturated candidates (with a message) instead
  of scoring them, and diagnoses an empty grid.
* `adl_joint_test()` matches the columns of a named `R` to the coefficients
  by name (they were applied by position). Wald tests stop with a clear
  message when the restrictions' covariance is singular.
* `reset_test()` rejects repeated powers and designs where the powers of the
  fitted values are collinear with the regressors.
* `comfac_test()`: the Wald form uses the term map and no longer runs the
  restricted optimizer; the LR optimizer checks convergence and no longer
  changes the caller's random-number state.
* `as_adl()` stops on `dynlm` options it cannot reproduce (`weights`,
  `subset`, ...) instead of silently refitting by OLS, treats bare regressors
  as lag 0, and checks that the converted fit reproduces the original.
* Documentation: the dynamic Chow forecast test is a plug-in test (no
  parameter uncertainty); Granger tests assume stationarity and
  "instantaneous causality" is a conditional contemporaneous association;
  long-run multipliers do not establish cointegration; the README states
  the package's scope.
* Requires R >= 4.2 (was 4.1): the suggested packages used by the tests
  (`dynlm` through `car` and `pbkrtest`) no longer install on R 4.1. CI
  checks R 4.2, the oldest supported version.

Fixes from the external audit of 2026-11-16 (high priority):

* Column names: a variable whose name equals a generated lag (e.g. `x_lag1`
  together with lag 1 of `x`) silently shared one column. Such clashes are
  now an error.
* Gaps in the time index (a missing period in quarterly dates, a skipped
  integer, ...) are now an error: lags are taken period by period, and a gap
  made the "lag" two periods old. Irregular data are not supported.
* `bg_test()` without an intercept could return a negative LM statistic; it
  now uses the uncentred R-squared and matches `lmtest::bgtest()`.
* `wu_hausman_test()` uses the fitted sample (never more observations than
  the fit; fewer only where the instrument lags need earlier data, with a
  message), resolves the regressor's column from the term map, checks the
  instrument rank, and reports the conventional first-stage F (and the HAC
  Wald / q with HAC); 10 is described as a rule of thumb.
* `weak_exogeneity_test()` stops on unstable models, conditions the marginal
  equation on the controls, the trend and (with several regressors) the
  other regressors, and documents its conditioning set and assumptions.
* Unstable models: `adl_long_run()` now returns `NA`, like `adl_dynamics()`,
  which still traces the multipliers over 20 lags by default.
* `adl_dynamics()`: the automatic horizon covers every distributed lag and
  stops only when the cumulative response has settled for a full AR cycle
  (it stopped at an early overshoot); the mean lag is documented as a signed
  moment.
* Lag values are validated before coercion (`lags = list(x = 1.9)` was
  silently truncated to 1).

* **Changed:** `start` now means the first observation of the estimation
  sample, with the lags taken from the observations before it, as in
  `dynlm` and EViews. In 0.1.0 and 0.2.0 the first max(lags) observations
  after `start` were lost to the lags.
* Exogeneity (ADL3): `wu_hausman_test()` (control-function test of
  predeterminedness, instruments from further lags, weak-instrument warning)
  and `weak_exogeneity_test()` (kappa = 0 in the regression for the change
  in x; `strong = TRUE` adds the Granger part). The regression includes
  x_{t-1} by default: without it the test is badly oversized for stationary
  variables (`x_level = FALSE` gives the ADL3 regression as written).
* `as_adl()` converts an existing `dynlm` fit into an `adl` fit.

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
