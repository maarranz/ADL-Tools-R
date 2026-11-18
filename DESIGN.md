# adltools — design notes

Decisions taken while planning the package (November 2026), revised after the
ADL1–ADL4 chapter revisions (14 November 2026). This file is excluded from the
built package (`.Rbuildignore`) and records the reasons behind the API, so
that later changes are deliberate.

## 1. Place in the ATSE suite

| Package | Role |
|:--|:--|
| arimatools | ARIMA estimation, diagnostics, outliers, dynamic forecasting |
| forecastingtools | Static forecasts, accuracy, forecast evaluation and comparison |
| unitroottools | Unit-root and stationarity tests |
| **adltools** | Stationary ADL/ECM models: specification, estimation, inference, causality, diagnostics, exogeneity, dynamic effects, dynamic forecasting |
| cointegration package (planned) | Everything about cointegration |

- adltools is the ADL counterpart of arimatools.
- Static (one-step, fixed-parameter) forecasts of ADL models and all forecast
  evaluation and comparison live in **forecastingtools**. That needs a
  `static_forecast()` method for ADL/lm fits and a
  `scheme = c("fixed", "recursive", "rolling")` argument there; it does not
  depend on adltools.
- Unit-root testing stays in **unitroottools**; adltools may depend on it
  (via `Remotes:` while it is GitHub-only).
- **Cointegration stays out**: bounds tests, ECM-based cointegration tests,
  FMOLS/DOLS, Gregory–Hansen, system (VECM) weak-exogeneity tests. The ECM
  *reparametrization* of a stationary ADL and the single-equation
  weak-exogeneity regression of ADL3 are in scope.
- Dependency direction: the cointegration package may depend on adltools and
  unitroottools, never the reverse.

## 2. Notation (as in ADL1–ADL4)

Labels in output follow the chapters:

| Symbol | Meaning |
|:--|:--|
| α₀, α_i | intercept and coefficients on y_{t−i} |
| β_j (one regressor), β_{ij} (regressor i) | coefficients on x_{t−j}, x_{i,t−j} |
| C(L), B(L), B_i(L) | lag polynomials; D_i(L) = B_i(L)/C(L) |
| k₀, k₁ (k_i) | long-run intercept and long-run multiplier(s), k_i = B_i(1)/C(1) |
| d_j | dynamic multipliers; M_q cumulative multipliers |
| b | speed of adjustment in the ECM, b = α₁ − 1 |
| λ, κ | coefficients of the marginal model for Δx; κ on the equilibrium error |
| ξ | Wu–Hausman control-function coefficient |
| ε_t | white-noise error of the ADL; u_t = ε_t / C(L) |

## 3. Model specification

### 3.1 Several regressors

Models may contain **n ≥ 1 explanatory variables**, ADL(p, q₁, …, q_n), each
with its own lag polynomial B_i(L), its own transfer function, dynamic
multipliers and long-run multiplier k_i. C(L) is common and alone governs
stability. Every function that is defined per regressor (long-run
multipliers, dynamic effects, exogeneity tests, COMFAC) handles all
regressors, reporting results per regressor where appropriate.

### 3.2 Internal representation: lag sets

The model is stored as a **named list of lag sets**, one per variable, e.g.

```r
lags = list(y = 1:2, x1 = 0:2, x2 = 1:3)   # x1 with x_t, x2 without (dead start)
```

Gaps are allowed (`y = c(1, 4)`). Missing lags are zeros in the
corresponding polynomial, so multipliers and forecasts work unchanged.

### 3.3 Three ways to specify a model, one representation

1. **Lag orders** (convenience):
   `adl(y ~ x1 + x2, data, order = c(p, q1, q2), first_lag = c(x1 = 0, x2 = 1))`.
   `order` gives maximum lags; `first_lag` the first lag of each regressor
   (default 0; `first_lag = 1` removes the contemporaneous term).
   `start` and `end` keep their usual time-series meaning (estimation
   sample), as in `ts()`, `window()` and `dynlm()`.
2. **Explicit lag sets**: `adl(y ~ x1 + x2, data, lags = list(...))`.
3. **dynlm-style formulas**: `y ~ L(y, 1:2) + L(x1, 0:2) + L(x2, 1:3)`, and
   existing `dynlm` fits.

All map to the lag-set list. Internal routines (Granger, Gets, lag selection,
COMFAC) loop over lag sets and never parse formulas. Each printed model shows
the equivalent dynlm formula.

### 3.4 Other terms

- **Unlagged regressors** after a bar: `y ~ x1 + x2 | X`. Accepts the
  intervention matrices of `arimatools::outlier_regressors()` (columns
  `AO2020Q2`, `LS2008M09`, `TC1990`), continued for forecasting with
  `arimatools::extend_interventions()`.
- **Deterministic terms**: intercept by default; `trend = TRUE` optional.
  **No seasonal dummies** for now.
- **Transformations** such as `log(y) ~ log(x)` are applied before lagging.

### 3.5 Structure detection

When a formula is a levels ADL (only `L(y, ·)`, `L(x_i, ·)`, intercept, trend,
unlagged regressors), the lag sets are parsed. Structural functions (dynamic
effects, long-run multipliers, dynamic forecasts) require that structure and
stop with a clear message for other forms (`d()`, ECM). Tests that only need
a fitted regression work with any formula.

### 3.6 The fitted object

The regression is built by the package (own lag matrix) on a common sample.
The result is an `lm` subclass recording the lag sets, so lmtest and sandwich
keep working and coefficients map exactly to α/β. Coefficient names are
readable: `y_lag1`, `x1`, `x1_lag1`.

### 3.7 Contemporaneous terms

Models with and without x_t must be equally easy to write (`first_lag`, lag sets,
formulas). Search procedures take a per-regressor rule
`contemporaneous = c("both", "include", "exclude")`:

- `include`: x_{i,t} always in the model;
- `exclude`: never (dead-start, Granger-type models);
- `both`: the search is run with and without x_{i,t} on the same sample and
  both results are reported side by side (with Wu–Hausman once available).

x_t is never dropped by a significance test alone: with simultaneity β̂₀ can
be significant because of the bias. The choice is the user's, as in ADL3.

## 4. Inference

- `adl_coefficients(fit)`: estimates, standard errors, t, p-values and
  confidence intervals, with α/β labels.
- `adl_joint_test(fit, terms)`: Wald test of linear restrictions, F and χ²
  (taxonomy restrictions, Gets steps).
- `adl_long_run(fit)`: k₀ and k_i for every regressor, with delta-method
  standard errors.
- `adl_forecast(fit, h, newdata, origin, level)` (v0.6.0): iterated
  forecasts row by row of the design (y lags replaced by forecasts; x, the
  unlagged terms and the trend from the data or `newdata`). Var e_h =
  σ² Σ_{i<h} ψ_i², ψ from 1/C(L); no parameter or x uncertainty (Lütkepohl
  §3.5 later). `newdata` dated (same index kind; rows after the data) or
  undated (rows follow the data; `next_index()` extends ts, yearqtr,
  yearmon, numeric and monthly/quarterly/daily dates). Origin inside the
  sample = in-sample dynamic simulation, flagged. Fan chart in base
  graphics, bands joined to the last observation.
- `j_test(fit_a, fit_b)` (v0.5.0): both directions on the common sample
  (dates matched), same dependent variable checked by value; `adl` or
  `dynlm` fits (a pure AR(1) is not an `adl`); refuses nested pairs (rank of
  [X_A, X_B] equal to that of one of them) and points to
  `adl_joint_test()`. z (default) or t; nonrobust or HAC. Outcome
  `neither` / `against_a` / `against_b` / `both` with the ADL3 wording.
  Checked against the two regressions written out.
- `adl_ecm(fit)` (v0.5.0): the ECM is linear in the ADL coefficients,
  ECM = Aθ + c (c = −1 for b), so its covariance is A V Aᵀ, exact and equal
  to that of the ECM regression (checked with nonrobust and HAC). b = −C(1),
  β*_i = B_i(1), γ_j = −Σ_{k>j} α_k, δ_i0 = β_i0, δ_ij = −Σ_{k>j} β_ik;
  absent lags are zero restrictions. Equilibrium-error form: b and k_i
  from `adl_long_run()`; the equilibrium error y − Σ k_i x_i as in
  `weak_exogeneity_test()`.
- `adl_dynamics(fit)` (v0.5.0): delta-method standard errors for d_h, M_h,
  k_i and the mean lag; the Jacobian of d_h follows the multiplier
  recursion (∂d_h/∂α_k = d_{h−k} + Σ α_i ∂d_{h−i}/∂α_k). The median lag is
  an integer (first crossing of one half), so its interval is the
  α/2, 1 − α/2 quantiles of the median lag over `nsim` draws from
  N(θ̂, V̂), with the share of draws where it is undefined. Checked against
  a numerical Jacobian and by Monte Carlo (ADL(1,1), T = 200: coverage
  0.95 for d_h, M_h, k; 0.93 for the mean lag; conservative for the median
  lag). `se = FALSE` gives the point estimates only.
- All take `cov_type = c("nonrobust", "HAC")` and `maxlags`.
- Asymptotic distributions by default (normal, χ²), with either
  covariance: the regressions are dynamic, so exact t and F do not hold.
  `dist = "t"` / `"F"` gives the EViews/`summary.lm()` figures. Differs from
  forecastingtools, where nonrobust means t/F.

## 5. Granger non-causality (ADL3)

Follows the direct (Sargent, 1976) procedure of ADL3, not
`lmtest::grangertest()`.

- **Four models**, each with its own lag orders (own lags, other variable's
  lags):
  1. y on lags of y and x, **no x_t** — tests x → y;
  2. x on lags of x and y, **no y_t** — tests y → x;
  3. model 1 plus x_t — instantaneous causality (y equation);
  4. model 2 plus y_t — instantaneous causality (x equation).
- **Two functions, one engine**:
  - `granger_test(y, x, p)`: the same order p in all four models.
  - `granger_test_general(y, x, ...)`: separate orders per model, supplied by
    the user or selected automatically (`select`, `max_order`: BIC, AIC, HQ
    over the (s, r) grid, or Gets; done in 0.4.0). The
    instantaneous models default to the orders of the corresponding Granger
    models. Unit test: with all orders equal to p it reproduces
    `granger_test()` exactly.
- **Common sample** for all models.
- **Distributions**: both χ² and F reported; `dist = c("chisq", "F")`
  decides the classification only. χ² is the default: the regressions are
  dynamic, so the exact F does not hold even with nonrobust errors.
- **Covariance**: `cov_type = "nonrobust"` by default (matches EViews and
  lmtest; the s = r identity holds), HAC optional. The identity is an OLS
  result; with HAC the two instantaneous statistics differ (stated in ADL3).
- **Classification**: one of the eight causal types, with symbol and wording
  from ADL3; flagged as ambiguous when the two instantaneous tests disagree.
- **Conditioning** (v0.6.0): `y ~ x | z1 + z2`; lags 1..`cond_lags`
  (default p, or the largest order) of every z in all four models, no z_t;
  output as before, with the classification labelled as conditional on
  the z's. Not with `select` (the IC grid would also range over the z
  lags).
- Stationarity assumed; stated in the documentation and one printed line.
- Validation: F equal to `lmtest::grangertest()` on the same sample;
  simulated bivariate VARs for each causal type.

## 6. COMFAC (ADL2), generalized

- Order-r common factor: C(L) and every B_i(L) share a degree-r factor ω(L).
  The model then reduces to an ADL(p − r, q₁ − r, …, q_n − r) with AR(r)
  errors: φ(L)y_t = c + Σ δ_i(L) x_{it} + u_t, ω(L)u_t = ε_t.
- **r · n restrictions** (n regressors); requires p ≥ r and q_i ≥ r.
- `comfac_test(fit, r = 1, method = c("LR", "Wald"))`:
  - **LR (default)**: unrestricted ADL by OLS against the restricted model
    with AR(r) errors, fitted by least squares on the same design matrix
    (same sample, no start-up difference); LR = T log(RSS_r / RSS_u) ~
    χ²(r·n). Given ω the restricted model is linear in (c, φ, δ) after
    filtering y and x by ω(L), so the RSS is concentrated on ω and
    minimized by a grid over the stationary region plus local refinement
    (checked against full NLS from random starts). Invariant, handles
    complex roots.
  - **Wald** (general since v0.7.0): Sargan's (1980) determinantal form.
    Polynomials ψ_j = C and the B_i without leading zeros, ordered by
    degree (C first among equals). Type 1 (Σ(f_j − r) > (n−1)(f−r)):
    rows of Ψ(k), eq. (8), first f + k + n − r of them; Type 2: Ψ_b (rows
    of Ψ₀ by decreasing right zeros), surplus vanishing rows dropped.
    Elimination row by row, tracking each reduced row's support (a pivot
    brings in its columns); a row with r elements left is vanishing; each
    element gives one determinant (pivot rows + row, pivot columns + its
    column) = algorithmic value × pivots. r·n determinants; derivatives
    exact (cofactors). ADL(1,1): α₁β₀ + β₁ (same statistic as the v0.2
    formula); one regressor, r = 1: the Sylvester resultant.
  - Checks: 63 random specifications (25 Type 1, 38 Type 2; 1–3
    regressors, p ≤ 4, dead-start regressors, r ≤ 3): constraints zero at
    an exact common factor, Jacobian of full rank r·n. Size (5% level,
    static regression with AR errors): ADL(2,2), r = 1: 0.032 (T = 200,
    2000 reps), 0.051 (T = 1000), 0.048 (T = 5000); ADL(3,3), r = 2: 0.030,
    0.042; ADL(2,2,2), r = 1: 0.030, 0.056. LR at T = 200: 0.04–0.07.
  - Not yet checked against PcGive output (Sargan's own example uses data
    not available to us). Bordering the leading square alone fails with
    two or more regressors (the square is singular under H0), which is why
    the vanishing-row algorithm is needed.
  - Output includes the restricted model's ω̂ (AR error coefficients).
  - Sequential testing r = 1, 2, … as in Sargan (1980).
- Documentation notes that Wald tests of nonlinear restrictions are not
  invariant to how the restriction is written (Gregory and Veall, 1985).

## 7. Diagnostics and stability (ADL3, ADL4)

- Breusch–Godfrey, ARCH-LM (no HAC: tests *for* the problem), RESET,
  Jarque–Bera.
- Chow breakpoint test.
- `chow_forecast_test(fit, T1, type = c("static", "dynamic"))`, parameters
  estimated once on T₁:
  - `static`: one-step forecasts with actual lagged values; exact
    Ω̂ = σ̂²(I + X₂(X₁′X₁)⁻¹X₂′) (EViews-compatible); σ̂²I gives the
    asymptotic version.
  - `dynamic`: a single chained path from T₁ (Hendry's predictive-failure
    test), Ω̂ built from the ψ-weights of 1/C(L).
- Hansen (1992): per coefficient, variance and joint (strucchange gives the
  joint statistic only, as far as known — check).

## 8. Exogeneity (ADL3)

For each regressor x_i:

- `wu_hausman_test(fit, instruments = NULL)`: control-function version
  (predeterminedness); default instruments x_{i,t−2}, y_{t−2}; tests ξ = 0.
- `weak_exogeneity_test(fit, lags = )`: builds ec_{t−1} = y_{t−1} − Σ k̂_i x_{i,t−1},
  regresses Δx_{i,t} on lags of Δx and Δy plus ec_{t−1}, tests κ = 0. Also
  reports the **strong-exogeneity** joint test (κ and all lagged Δy terms).
- `super_exogeneity_test()`: breaks in the marginal model of x_i summarized
  by impulse/step dummies (from the outlier/break tools), then those dummies
  tested in the conditional model (Engle and Hendry, 1993).

## 9. Outliers (v0.3.x)

Decided 2026-11-17:

- `adl_outliers(fit, types, cval, method)`: impulse and step dummies **in the
  ADL equation** (as `gets::isat()` and `super_exogeneity_test()`), not
  outlier types defined on y. An equation impulse is an innovation outlier;
  a step shifts the intercept, so y moves gradually to a new level. An AO in
  y shows as an impulse at τ and a smaller opposite one at τ+1; an immediate
  LS in y as a step plus an impulse. Documented in the help page.
- Engines: own forward-backward search on t-statistics (default; Chen-Liu
  type adapted to OLS, transparent for teaching) and `gets::isat()` (IIS/SIS)
  as `method = "isat"`, `gets` in Suggests.
- Critical value: the tsoutliers rule (3 for T <= 50, 4 for T >= 450,
  3 + 0.0025(T - 50) in between); isat uses the matching two-sided p-value.
- Output: the refitted ADL (dummies as unlagged regressors, named
  `I_2005Q1`, `S_2005Q1`), and the dates as lists for
  `super_exogeneity_test(impulses =, steps =)`.
- Not done: arimatools-format names (AO/LS), lagged intervention columns,
  outliers in x as a separate type, pure AR models (adl() needs a regressor;
  use `adl(x ~ y, first_lag = 1)` for the marginal model).

## 10. Automatic modelling (Gets)

- **Tier 1, `adl_gets()`** (done in 0.3.0.9000) — the ADL3 procedure:
  general ADL(p, …, p) on the common sample (p by frequency: 4, 12, else 2);
  diagnostics on the general model (warning if it fails); drop the highest
  lags jointly while not rejected; stop at the first rejection; back up if a
  simplified model fails a diagnostic it did not fail before (BG of order p,
  ARCH, RESET by default; JB optional); report the full path. Uses the
  `contemporaneous` rule (§3.7), with `"include"` as default. No
  `method = "ic"`: information criteria stay in `adl_ic()`. Monte Carlo
  (true ADL(2,2), T = 200, 300 reps): selects 2 in 82%, 3 or 4 in 18%,
  never fewer.
  Shares its engine with automatic lag selection in
  `granger_test_general()` and the sequential COMFAC search.
- **Tier 2, `adl_getsm(fit)`** (v0.7.0) — the GUM's design matrix to
  `gets::arx()` (mconst for the intercept), `gets::getsm()` with `keep` =
  intercept, unlagged terms, trend and (default) every x_t; t.pval =
  `alpha`; Ljung-Box AR/ARCH diagnostics at 2.5% (gets' default); terminals
  by BIC (`ic`). Selected columns refitted with `adl(lags = ...)` on the
  GUM's data and `max_lag` (same sample; gaps allowed); coefficients equal
  getsm's. If every regressor is dropped, `fit = NULL`. `gets` in
  Suggests.
- Documentation notes that post-selection standard errors are conditional on
  the search.

## 11. Conventions shared with the suite

- **Names**: snake_case; tests end in `_test`. The `adl_` prefix only where a
  name would mask an arimatools function (`adl_dynamic_effects`,
  `adl_dynamic_forecast`, `plot_adl_inverse_roots`).
- **Test results**: an `adltools_test` list with the fields of
  forecastingtools (`statistic`, `p.value`, `distribution`, `df`/`df_num`/
  `df_denom`, `null_hypothesis`, `reject_null`, `alpha`, `cov_type`,
  `maxlags`, `nobs`, `coefficients`) and the same printed layout, always
  with a `Covariance:` line.
- **HAC**: computed as in forecastingtools, `sandwich::NeweyWest(fit, lag = L,
  prewhite = FALSE, adjust = FALSE)`, default bandwidth
  min(n − 1, ⌊4(n/100)^(2/9)⌋). Default in adltools: `"nonrobust"`.
- **Tables**: data frames. **Plots**: ggplot2/patchwork, `interactive = TRUE`
  for plotly.
- **Data types**: numeric vectors, `ts`, `zoo`/`xts`, data frames/tibbles with
  a time column, tsibbles — following unitroottools, extended to several
  series on a common time index. tibble, tsibble, zoo, xts in Suggests;
  `tidy()`/`glance()` via generics.
- **Tests**: every function checked against a known result (EViews output,
  lmtest/strucchange, or a simulated DGP).

## 12. Data handling

- **Missing values**: leading and trailing NAs are trimmed with a message;
  internal NAs stop with an error (dropping them would break the lag
  structure).
- **Time index for data frames/tibbles**: a time column given by `time =`, or
  detected when exactly one Date/yearmon/yearqtr column exists; otherwise row
  order, with a message.
- **Sample**: `start`/`end` select the estimation sample, as in `dynlm` and
  EViews: `start` is the first observation estimated and the lags use the
  observations before it (decided 2026-11-16; 0.1.0 and 0.2.0 used the
  first max(lags) observations after `start` instead). Models with
  different lags and the same `start` share one sample. Without `start`, the
  first max(lags) (or `max_lag`) observations of the data are the initial
  values; if there are too few observations before `start`, the sample
  starts as early as possible, with a message.

## 13. Implementation plan

**Dependencies.** Imports: stats, zoo (index alignment, lagging), sandwich
(HAC), lmtest (`bgtest()`, `resettest()`), generics (`tidy`, `glance`).
Suggests: testthat, tibble, tsibble, xts, dynlm and ARDL (cross-checks in
tests), knitr, rmarkdown.

**v0.1.0 files**

| File | Contents |
|:--|:--|
| `R/input.R` | `as_adl_data()` (internal): any input type to a zoo matrix with regular index; NA handling; names |
| `R/lags.R` | `order`/`first_lag`, `lags` or formula to the lag-set list; `build_design()` (design matrix, readable names, common sample) |
| `R/formula.R` | Parser for `y ~ x1 + x2 \| D` and `y ~ L(y, 1:2) + L(x, 0:2) \| D`; equivalent dynlm formula for printing |
| `R/adl.R` | `adl()`, `print`, `summary`; object stores lag sets, sample, index, data |
| `R/vcov.R` | `adl_vcov(fit, cov_type, maxlags)`, identical to forecastingtools |
| `R/test-class.R` | `adltools_test` constructor and print method |
| `R/inference.R` | `adl_coefficients()`, `adl_joint_test()`, `adl_long_run()` |
| `R/granger.R` | engine, `granger_test()`, `granger_test_general()` (`select =` added in v0.4.0), classification |
| `R/tidy.R` | `tidy()`/`glance()` methods |
| `vignettes/granger.Rmd` | Granger causality following ADL3 |
| `data/` | small simulated ADL(1,1) dataset with known parameters |

Order of work (one commit each): input and lags → `adl()` → vcov and test
class → inference → Granger → tidy/glance → vignette and data → release
(tag `v0.1.0`).

**v0.2.0 files**

| File | Contents |
|:--|:--|
| `R/diagnostics.R` | `bg_test()`, `reset_test()` (lmtest underneath); `arch_test()`, `jb_test()` written directly (no FinTS/tseries); `adl_diagnostics()` in ADL3's order |
| `R/chow.R` | `chow_break_test()` (dummy interactions, HAC possible); `chow_forecast_test()` static (dummy-variable regression = exact Ω̂ form) and dynamic (Ω̂ from ψ-weights of 1/C(L), general p) |
| `R/comfac.R` | `comfac_test(fit, r, method)`: restricted model parametrized by (c, φ, δ_i, ω), mapped to unrestricted coefficients by polynomial multiplication, RSS minimized with `optim()` on the same design matrix (same sample by construction); start values from the roots of Ĉ(L) |

Outlier detection: `adl_outliers()` (section 9).

**v0.5.0 files**

| File | Contents |
|:--|:--|
| `R/dynamics.R` | `adl_dynamics()` with delta-method standard errors (Jacobian by the multiplier recursion) and the simulated median-lag interval |
| `R/ecm.R` | `adl_ecm()`: ECM = Aθ + c, exact covariance A V Aᵀ; equilibrium-error form and the equilibrium error |
| `R/jtest.R` | `j_test()`: both directions on the common sample; `adl` or `dynlm` fits; nesting check |

**v0.6.0 files**

| File | Contents |
|:--|:--|
| `R/forecast.R` | `adl_forecast()`, `next_index()` (extends a regular index), print and fan-chart plot |
| `R/granger.R` | formula interface and conditioning variables (`parse_granger_formula()`, `cond_lags`) |

**v0.7.0 files**

| File | Contents |
|:--|:--|
| `R/comfac.R` | `sargan_structure()` (rows of Ψ(k) or Ψ_b, elimination with vanishing rows, the determinants), `sargan_constraints()` (values and exact Jacobian) |
| `R/getsm.R` | `adl_getsm()`: adapter to `gets::arx()`/`getsm()`, refit as `adl` |

**Testing.** Development tests run in a Linux R installation in the
assistant's workspace; final `R CMD check --as-cran` on the M1 Pro and CI.

## 14. Roadmap

- **v0.1.0** (released 2026-11-16) — specification and estimation (`adl()` from `order`/`first_lag`,
  explicit `lags` or its own dynlm-style formulas; several regressors;
  unlagged dummies; base and tidy inputs); inference (`adl_coefficients()`,
  `adl_joint_test()`, `adl_long_run()`); Granger (`granger_test()`,
  `granger_test_general()` with user-supplied orders only); `tidy()`/
  `glance()`; dynamic effects (`adl_dynamics()`: multipliers, mean and
  median lag, inverse roots; point estimates); `adl_ic()` table of
  information criteria on a common sample; one vignette on Granger
  causality.
- **v0.2.0** (released 2026-11-16) — diagnostics (BG, ARCH, RESET, JB), Chow breakpoint,
  `chow_forecast_test()` (static and dynamic), `comfac_test()` (LR, and Wald
  for the ADL(1,1)).
- **v0.2.1** (released 2026-11-16) — `as_adl()` for existing `dynlm` fits;
  `start`/`end` as in dynlm/EViews; exogeneity (`wu_hausman_test()`,
  `weak_exogeneity_test()` with strong exogeneity); fixes from the first
  audit; R >= 4.2.
- **v0.2.2** (released 2026-11-17) — workflow vignette; fixes from the
  second audit.
- **v0.3.0** (released 2026-11-17) — stability without a break date
  (`recursive_residuals()`, `cusum_test()`, `cusumsq_test()` with plots,
  `hansen_test()` per parameter and joint); `super_exogeneity_test()`.
- **v0.4.0** (released 2026-11-17) — `adl_outliers()` (own search or
  `gets::isat()`), `adl_gets()` tier 1, automatic order selection in
  `granger_test_general()` (`select = "bic"/"aic"/"hq"/"gets"`), the
  `slopes` option of `super_exogeneity_test()`; fixes from the third audit.
- **v0.5.0** (released 2026-11-18) — standard errors in `adl_dynamics()`
  (delta method; simulated interval for the median lag); `adl_ecm()`;
  `j_test()`.
- **v0.6.0** (released 2026-11-18) — `adl_forecast()` (dynamic forecasts
  on a known path of x, ψ-weight variance, fan charts); Granger with
  conditioning variables (`y ~ x | z`) and a formula interface.
- **v0.7.0** (released 2026-11-18) — Sargan's (1980) determinantal COMFAC
  Wald test for any r and n; `adl_getsm()` (tier 2 Gets through the `gets`
  package).
- **Later** — parameter uncertainty in forecast variances (Lütkepohl
  2005, section 3.5).
- **Elsewhere** — forecastingtools: `static_forecast()` for ADL fits and the
  `scheme` argument.

## 15. Open questions

1. Test data: simulated DGPs for unit tests; which real dataset to
   cross-check against EViews?
2. ~~ECM view~~: decided, `adl_ecm(fit)` (v0.5.0).
3. ~~Outlier critical values~~: decided, the tsoutliers rule (section 9).
4. ~~Granger with several regressors~~: decided, `y ~ x | z` (v0.6.0).

## 16. References (checked against publisher or RePEc records, Nov 2026)

- Engle, R. F., Hendry, D. F. and Richard, J.-F. (1983). Exogeneity.
  *Econometrica* 51(2), 277–304.
- Engle, R. F. and Hendry, D. F. (1993). Testing superexogeneity and
  invariance in regression models. *Journal of Econometrics* 56 (pages 119–139
  not verified; check against `referencesADL.bib`).
- Gregory, A. W. and Veall, M. R. (1985). Formulating Wald tests of nonlinear
  restrictions. *Econometrica* 53(6), 1465–1468.
- Hendry, D. F. and Mizon, G. E. (1978). Serial correlation as a convenient
  simplification, not a nuisance: A comment on a study of the demand for
  money by the Bank of England. *Economic Journal* 88(351), 549–563.
- Mizon, G. E. and Hendry, D. F. (1980). An empirical application and Monte
  Carlo analysis of tests of dynamic specification. *Review of Economic
  Studies* 47(1), 21–45.
- Sargan, J. D. (1980). Some tests of dynamic specification for a single
  equation. *Econometrica* 48(4), 879–897.
- Sargent, T. J. (1976). A classical macroeconometric model for the United
  States. *Journal of Political Economy* 84(2), 207–237.

Other chapter references (Hansen 1992, Hendry 1979, Urbain 1992, Wu 1973,
Hausman 1978, Chow 1960, Davidson–MacKinnon 1981) are in `referencesADL.bib`.
