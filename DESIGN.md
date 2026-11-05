# adltools — design notes

Decisions taken while planning the package (November 2026). This file is
excluded from the built package (`.Rbuildignore`) and records the reasons
behind the API, so that later changes are deliberate.

## 1. Place in the ATSE suite

| Package | Role |
|:--|:--|
| arimatools | ARIMA estimation, diagnostics, outliers, dynamic forecasting |
| forecastingtools | Static forecasts, accuracy, forecast evaluation and comparison |
| unitroottools | Unit-root and stationarity tests |
| **adltools** | Stationary ADL/ECM models: estimation, inference, causality, diagnostics, dynamic effects, dynamic forecasting |
| cointegration package (planned) | Everything about cointegration |

- **adltools is the ADL counterpart of arimatools.**
- **Static (one-step) forecasts of ADL models, and all forecast evaluation and
  comparison, live in forecastingtools.** That needs a `static_forecast()`
  method for ADL/lm fits there; it does not depend on adltools.
- **Unit-root testing stays in unitroottools.** adltools may depend on it
  (via `Remotes:` while it is GitHub-only).
- **Cointegration stays out**: bounds tests, ECM-based cointegration tests,
  FMOLS/DOLS, Gregory–Hansen, weak exogeneity with respect to cointegrating
  parameters (Urbain-type). The ECM *reparametrization* of a stationary ADL is
  in scope.
- Dependency direction: the cointegration package may depend on adltools and
  unitroottools, never the reverse.

## 2. Model formulation

Two entry points, one internal representation.

- **dynlm-style formulas** for users: `y ~ L(y, 1:2) + L(x, 0:2)`. Any lag
  pattern, transformations and ECM forms can be written. Plain `dynlm` fits
  estimated by students are also accepted.
- **Lag orders** for internal routines and for convenience:
  `adl(y ~ x, data, order = c(p, q))`, several regressors
  `order = c(p, q1, q2)`.
  - `order` gives the maximum lag of each variable, as in ADL(p, q).
  - `start` gives the first lag of each regressor (default 0). `start = 1`
    drops the contemporaneous term (dead-start model), e.g.
    `start = c(x1 = 0, x2 = 1)`.
- **Unlagged regressors** (dummies, outlier/intervention variables) after a
  bar: `y ~ x | X`. They accept the matrices built by
  `arimatools::outlier_regressors()` (columns `AO2020Q2`, `LS2008M09`,
  `TC1990`) and are continued for forecasting with
  `arimatools::extend_interventions()`.
- **Deterministic terms**: intercept by default, `trend = TRUE` optional.
  **No seasonal dummies** for now.
- **Transformations** such as `log(y) ~ log(x)` are applied before lagging.
- The regression is built by the package (own lag matrix), and the result is
  an `lm` subclass recording p, q and the starting lags, so that lmtest and
  sandwich keep working and coefficients map exactly to α/β.
- Coefficient names are readable (`y_lag1`, `x`, `x_lag1`).
- **Structure detection**: when a formula is a levels ADL (only `L(y, ·)`,
  `L(x, ·)`, intercept, trend, unlagged regressors), p, q and the starting
  lags are parsed. Structural functions (dynamic effects, long-run
  multiplier, dynamic forecasts) require that structure and stop with a clear
  message for other forms (`d()`, ECM). Tests that only need a fitted
  regression work with any formula.
- The equivalent dynlm formula is printed with each model.

Open: ECM view as a separate `adl_ecm(fit)` or a `form = "ecm"` argument
(leaning towards estimation in ADL form, ECM as a view).

## 3. Granger non-causality

Follows the direct (Sargent, 1976) procedure in the chapters (ADL3, DynReg3),
not `lmtest::grangertest()`.

- **Four models**, each with its own lag orders (own lags, other variable's
  lags):
  1. y on lags of y and x, **no x_t** — tests x → y;
  2. x on lags of x and y, **no y_t** — tests y → x;
  3. model 1 plus x_t — instantaneous causality (y equation);
  4. model 2 plus y_t — instantaneous causality (x equation).
- **Two functions, one engine**:
  - `granger_test(y, x, p)`: the same order p in all four models, as in other
    software.
  - `granger_test_general(y, x, ...)`: separate orders per model, either
    supplied by the user or selected automatically (maximum lag and a
    criterion: AIC/BIC or Gets). The instantaneous models default to the
    orders of the corresponding Granger models.
  - Unit test: the general function with all orders equal to p reproduces
    `granger_test()` exactly.
- **Common sample**: all models use the same observations (the first
  max(lag orders) are dropped). Required for the s = r identity to hold
  exactly.
- **Distributions**: every test reports both F and χ² (qF), each with its
  p-value. `dist = c("F", "chisq")` only decides the classification.
- **Covariance**: `cov_type = "nonrobust"` by default (results match EViews
  and lmtest, and the chapters' identity holds); `cov_type = "HAC"` with
  `maxlags` optional. With HAC the statistics are Wald tests (χ² = W,
  F = W/q), and the instantaneous statistics of models 3 and 4 differ even
  with equal orders: the identity is an OLS result.
- **Classification**: one of the eight causal types, with symbol and wording
  from ADL3. When the two instantaneous tests disagree (general orders, or
  HAC), the result is flagged as ambiguous and both are shown.
- **Scope**: bivariate only; third variables later.
- **Stationarity**: assumed; stated in the documentation and in one printed
  line (see unitroottools). No automatic unit-root test.
- Multiple testing across the three tests is documented, not adjusted.

Validation: simple-causality F equal to `lmtest::grangertest()` on the same
sample; simulated bivariate VARs for each causal type (correlated errors for
instantaneous causality) are classified correctly in long samples.

## 4. Parameter inference

- `adl_coefficients(fit)`: estimates, standard errors, t, p-values and
  confidence intervals, with α/β labels for levels ADLs.
- `adl_joint_test(fit, terms)`: Wald test of joint zero restrictions, F and
  χ² (Gets steps, taxonomy restrictions such as the static model).
- `adl_long_run(fit)`: long-run multiplier with delta-method standard error.
- All take `cov_type` and `maxlags`.

## 5. Conventions shared with the suite

- **Names**: snake_case; tests end in `_test` (`granger_test`,
  `chow_forecast_test`, `hansen_test`). The `adl_` prefix is used only where
  a name would mask an arimatools function (`adl_dynamic_effects`,
  `adl_dynamic_forecast`, `plot_adl_inverse_roots`).
- **Test results**: an `adltools_test` list with the fields of
  forecastingtools (`statistic`, `p.value`, `distribution`, `df`/`df_num`/
  `df_denom`, `null_hypothesis`, `reject_null`, `alpha`, `cov_type`,
  `maxlags`, `nobs`, `coefficients`) and the same printed layout, always
  including a `Covariance:` line.
- **HAC**: `cov_type = "HAC"` / `"nonrobust"` and `maxlags`, computed exactly
  as in forecastingtools: `sandwich::NeweyWest(fit, lag = L, prewhite = FALSE,
  adjust = FALSE)`, default bandwidth min(n − 1, ⌊4(n/100)^(2/9)⌋).
  Default in adltools: `"nonrobust"` (in forecastingtools HAC is the default
  because h-step errors overlap; here it is "yes, if in doubt", ADL3).
- **Tables**: data frames. **Plots**: ggplot2/patchwork, `interactive = TRUE`
  for plotly (works in both HTML and Typst output of the chapters).
- **Data types (base and tidy)**: input handling follows unitroottools —
  numeric vectors, `ts`, `zoo`/`xts`, data frames/tibbles with a time column,
  tsibbles — extended to several series on a common time index. tibble,
  tsibble, zoo and xts in Suggests. `tidy()` and `glance()` methods via the
  generics package.

## 6. Outliers (v0.2)

- Types: AO, LS, TC, IO and outliers in x. Documentation presents AO and LS
  first; IO and outliers in x as advanced material.
- In an ADL an unlagged dummy passes through 1/α(L): an AO dummy alone acts
  as an innovational outlier, an LS dummy gives a gradual shift. A genuine
  additive outlier in y needs the dummy and its lags up to p. Documented,
  with an option to lag intervention columns.
- Detection: Chen–Liu-type search adapted to OLS (dummy signatures per type),
  re-tested jointly. Compare with `gets::isat()` (IIS/SIS) before writing.
- An IO in an ADL(1, q) behaves like a TC with δ = α₁: when the estimate of
  α₁ is close to δ the two are nearly indistinguishable, and the output says
  which competing type came close.
- Output in the arimatools format, so detected outliers go straight back into
  `y ~ x | X` and `extend_interventions()`.
- Open: fixed critical values (3–3.5) or IIS-style levels tied to T.

## 7. Roadmap

- **v0.1** — formulation, estimation, inference, causality:
  `adl()` and parsing of dynlm fits; `adl_coefficients()`,
  `adl_joint_test()`, `adl_long_run()`; `granger_test()`,
  `granger_test_general()`; base and tidy inputs.
- **v0.2** — diagnostics (Breusch–Godfrey, ARCH-LM, RESET, Jarque–Bera,
  Chow breakpoint, Chow forecast) and outlier detection.
- **Later** — Hansen (1992) parameter instability (per coefficient, variance
  and joint; strucchange gives the joint statistic only), Durbin–Wu–Hausman
  (control function), J-test with the four-outcome reading, Gets and
  automatic lag selection, dynamic effects and long-run multiplier with
  bands, inverse roots, ECM view, dynamic forecasting with fan charts.
- **Elsewhere** — forecastingtools: `static_forecast()` for ADL fits.

## 8. Open questions

1. ECM: `adl_ecm(fit)` or a `form` argument.
2. Test data for examples and unit tests (chapter datasets, simulated DGPs,
   a real dataset cross-checked against EViews).
3. Outlier critical values; whether to detect outliers in x by default.
