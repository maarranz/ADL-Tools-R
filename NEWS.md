# adltools 0.7.1

Fixes from the fourth external audit (2026-11-18).

* `adl_ecm()`: a term is shown as fixed (zero standard error, no test) only
  when it depends on no coefficient (b = -1 without lags of y). Before,
  any standard error below a threshold was set to zero, so rescaling x
  (e.g. by 1e10) gave point intervals. Row labels are now unique when a
  control's name coincides with a generated ECM term (e.g. `x_lag1`).
* `adl_forecast()`: rows of `newdata` after the end of the data are
  selected with half the data's spacing, not an epoch-relative tolerance
  that dropped date-times one second apart. Dates extended for monthly
  data on a day capped at short months (the 30th, 28th in February) keep
  that convention. The fan chart works for any set of levels (column tags
  are built the same way for the table and the plot).
* Time indexes: regular sub-second grids are no longer rejected as
  irregular (the tolerance allows for the floating-point resolution of
  epoch seconds); gap messages show milliseconds.
* `adl_getsm()` and `adl_outliers()`: the refitted model keeps
  non-syntactic variable names (formulas are rebuilt with quoted names).
* `j_test()`: rank-deficient models (on the common sample) and `dynlm`
  fits with aliased coefficients, weights or an offset are refused with a
  clear message, instead of failing with "0 (non-NA) cases".
* `cusum_test()`: very small `alpha` (e.g. 1e-50) works; the p-value is
  computed on the log scale.
* `adl_outliers()`: samples too short for a search stop with a message,
  and a search without candidates returns no outliers.
* `adl_dynamics()`: `max_horizon` must be a single integer >= 0.
* Documentation: in `granger_test_general(select = )` the orders are
  chosen on the sample of the largest model, and the tests run on the
  common sample of the selected orders.

# adltools 0.7.0

* `adl_getsm()`: term-by-term, multi-path general-to-specific search
  through `gets::getsm()` (tier 2 of the Gets design). Takes a general ADL
  fit, searches its terms one by one (so a lag can stay without the ones
  before it), and refits the selection as an `adl` on the same sample.
  The intercept, unlagged terms, the trend and (by default) every x_t are
  never tested; `keep` adds others. `adl_gets()` remains the lag-by-lag
  procedure of ADL3.
* `comfac_test(method = "Wald")` is now Sargan's (1980) determinantal Wald
  test, for any order r and any number of regressors (it was limited to the
  ADL(1,1)). The common-factor hypothesis is written as r * n determinants
  of shifted coefficient rows, following Sargan's elimination algorithm
  (Type 1 and Type 2 cases, vanishing rows), with exact derivatives. For the
  ADL(1,1) it is the chapter's alpha_1 beta_0 + beta_1 = 0, with the same
  statistic as before. The help page explains Sargan's warning on reading
  the tests for successive r. The determinants are returned in
  `constraints`.

# adltools 0.6.0

* `granger_test()` and `granger_test_general()` accept a formula, `y ~ x`,
  and conditioning variables, `y ~ x | z1 + z2`: lags 1 to `cond_lags` of
  each z enter all four regressions (never z_t), so the tests and the
  classification are conditional on the past of the z's. The name-based
  calls are unchanged. Order selection with `select` is not available with
  conditioning variables. New section in `vignette("granger")`: a common
  cause that the conditioning removes.
* `adl_forecast()`: dynamic (iterated) forecasts of y from an ADL, as in
  ADL4, with the regressors on a known path: from the data inside the
  sample (with the actual values for comparison) and from `newdata` beyond
  it (dated, or rows following the data). Standard errors from the
  psi-weights of 1/C(L), sigma^2 sum psi_i^2, without parameter
  uncertainty; normal intervals at several levels and a fan chart
  (`plot()`). `tidy()` and `glance()` methods.

# adltools 0.5.0

* `j_test()`: the J-test of Davidson and MacKinnon (1981) for two
  non-nested models (ADL3), run in both directions on their common sample,
  with the four-outcome reading of the chapter. Takes `adl` fits or
  `dynlm` fits (so an AR(1) can be compared with a distributed lag), with
  conventional or HAC standard errors. Nested models are refused, with a
  pointer to `adl_joint_test()`. `tidy()` and `glance()` methods.
* `adl_ecm()`: the error correction form of a fitted ADL (ADL2), for any
  orders and several regressors: the adjustment coefficient b, the
  B_i(1) on the lagged levels and the coefficients of the lagged
  differences, with exact standard errors (the ECM is a linear
  reparametrization, so they equal those of the ECM regression, with
  either covariance). Also the equilibrium-error form (b and the long-run
  multipliers k_i) and the equilibrium error as a series. `tidy()` and
  `glance()` methods.
* `adl_dynamics()` reports standard errors: delta method for the dynamic and
  cumulative multipliers, the long-run multiplier and the mean lag (with
  `cov_type`, `maxlags`, `alpha` and `dist` as in the other inference
  functions), and a simulated interval for the median lag, which is an
  integer and has no derivative (`nsim` draws from the asymptotic
  distribution of the coefficients; the share of draws where it is
  undefined is reported). `plot()` draws the pointwise confidence band
  (`ci = FALSE` to omit it). `se = FALSE` gives the previous output.

# adltools 0.4.0

* `adl_outliers()`: impulse and step dummies in the ADL equation, found by a
  forward-backward search on t-statistics (default) or by indicator
  saturation with `gets::isat()` (`method = "isat"`), with the tsoutliers
  critical value (3 to 4, growing with T). Returns the ADL refitted with the
  dummies, and the dates in the form `super_exogeneity_test()` takes.
* `adl_gets()`: the general-to-specific procedure of ADL3. From a general
  ADL(p, p) on a common sample (p = 4 quarterly, 12 monthly by default), the
  highest lags of y and every x are dropped jointly while their Wald test
  does not reject; the diagnostics (BG, ARCH, RESET by default) are re-run
  on every simpler model, and the search backs up if one fails. x_t is never
  dropped by significance: `contemporaneous = "include"`, `"exclude"` or
  `"both"`. The full path is printed.
* `granger_test_general(select = )`: automatic order selection for each
  Granger equation, by BIC, AIC or HQ over all (s, r) up to `max_order`, or
  by the general-to-specific search (documented as a pre-test). The
  instantaneous models use the selected orders.

## Fixes from the third audit

* `super_exogeneity_test()`: a control whose name coincided with an internal
  dummy column (e.g. `dummy.1`) was tested in place of the break dummy. The
  internal names now avoid every column of both regressions. New argument
  `slopes = TRUE` also interacts the steps with the regressors, because
  additive dummies can miss a change in a slope.
* Times (`POSIXct`): a missing minute or hour was accepted as consecutive
  (the check used dates only). Steps are now checked in seconds, with
  calendar periods only when every observation has the same clock time
  (daily data across daylight-saving changes stay regular).
* Dates of dummies (`super_exogeneity_test()`, `adl_outliers()`) are matched
  to exactly one row; the tolerance was large enough to match several rows
  on a seconds grid.
* `cusumsq_test()`: the decision now follows the Monte Carlo p-value, and
  the bound is the order statistic that agrees with it; `nsim` must satisfy
  (nsim + 1) * alpha > 1. Before, a tiny `nsim` could reject with p = 0.5.

# adltools 0.3.0

* Stability tests (ADL3): `recursive_residuals()`, `cusum_test()` and
  `cusumsq_test()` (Brown, Durbin and Evans, 1975), with `plot()` methods
  drawing the bounds, and `hansen_test()` (Hansen, 1992), which tests each
  coefficient and the error variance, and all jointly, without a break date.
  CUSUM uses the EViews scaling (standard error of the regression); the
  CUSUM-of-squares bounds come from the exact distribution by simulation.
* `super_exogeneity_test()` (Engle and Hendry, 1993, as outlined in ADL3):
  impulse and step dummies for the breaks in the marginal process of x are
  tested in the regression for the change in x (they should matter), in the
  ADL (they should not: invariance) and with weak exogeneity, in one table.
  Warns when the marginal model shows no break, so the test is uninformative.

# adltools 0.2.2

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
