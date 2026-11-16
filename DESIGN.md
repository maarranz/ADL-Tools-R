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
    the user or selected automatically (maximum lag; AIC/BIC or Gets). The
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
- **Scope**: bivariate for now; conditioning on other regressors (third
  variables) later.
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
    with AR(r) errors by conditional least squares (`stats::arima(...,
    method = "CSS")`), same sample; LR = T log(RSS_r / RSS_u) ~ χ²(r·n).
    Invariant, handles complex roots.
  - **Wald**: for the ADL(1,1), the ADL2 formula α₁β₀ + β₁ = 0, so students
    can reproduce the chapter. General Wald later, in Sargan's determinantal
    form (no choice of root; comparable with PcGive).
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

## 9. Outliers (v0.2.x)

- Types: AO, LS, TC, IO, and outliers in x. Documentation presents AO and LS
  first; IO and outliers in x as advanced.
- In an ADL an unlagged dummy passes through 1/C(L): an AO dummy alone acts
  as an IO; an LS dummy gives a gradual shift; a genuine AO in y needs the
  dummy and its lags up to p. Option to lag intervention columns.
- An IO in an ADL(1, q) behaves like a TC with δ = α₁: the output says which
  competing type came close.
- Candidate engine: `gets::isat()` (impulse/step indicator saturation), to be
  checked; otherwise a Chen–Liu-type search adapted to OLS.
- Output in the arimatools format, so detected outliers go straight back into
  `y ~ x | X` and `extend_interventions()`.

## 10. Automatic modelling (Gets)

- **Tier 1, `adl_gets()`** — the ADL3 procedure: general ADL(p, …, p) on the
  common sample; diagnostics on the general model; drop the highest lags
  jointly while not rejected; stop at the first rejection; back up if a
  simplified model fails a diagnostic; report the full path. `method = "ic"`
  for AIC/BIC over the same grid. Uses the `contemporaneous` rule (§3.7).
  Shares its engine with automatic lag selection in
  `granger_test_general()` and the sequential COMFAC search.
- **Tier 2** — term-by-term multi-path search through an adapter to the
  `gets` package (`arx()`, `getsm()`), returning an `adl` fit. `gets` in
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
- **Sample**: `start`/`end` select the estimation sample; the first
  max(lags) observations of that sample are used as initial values, so all
  models built from the same data and lag bound share one common sample.

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
| `R/granger.R` | engine, `granger_test()`, `granger_test_general()` (`select =` stops with "not yet implemented"), classification |
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

Outlier detection moves to **v0.2.x** (depends on the `gets::isat()`
decision).

**Testing.** Development tests run in a Linux R installation in the
assistant's workspace; final `R CMD check --as-cran` on the M1 Pro and CI.

## 14. Roadmap

- **v0.1.0** — specification and estimation (`adl()` from `order`/`first_lag`,
  explicit `lags` or its own dynlm-style formulas; several regressors;
  unlagged dummies; base and tidy inputs); inference (`adl_coefficients()`,
  `adl_joint_test()`, `adl_long_run()`); Granger (`granger_test()`,
  `granger_test_general()` with user-supplied orders only); `tidy()`/
  `glance()`; one vignette on Granger causality.
- **v0.1.x** — accept existing `dynlm` fits.
- **v0.2.0** — diagnostics (BG, ARCH, RESET, JB), Chow breakpoint,
  `chow_forecast_test()` (static and dynamic), `comfac_test()` (LR, and Wald
  for the ADL(1,1)).
- **v0.2.x** — outlier detection.
- **v0.3** — exogeneity (`wu_hausman_test()`, `weak_exogeneity_test()` with
  strong exogeneity, `super_exogeneity_test()`), `adl_gets()` tier 1, and
  automatic lag selection in `granger_test_general()` (same engine; the
  argument exists from v0.1.0 and stops with "not yet implemented").
- **Later** — Hansen (1992) per coefficient; J-test with the four-outcome
  reading; dynamic multipliers, mean and median lag, inverse roots; ECM view;
  dynamic forecasts with ψ-weight variance and fan charts; Sargan's general
  COMFAC Wald; `gets` adapter (tier 2); Granger with conditioning variables.
- **Elsewhere** — forecastingtools: `static_forecast()` for ADL fits and the
  `scheme` argument.

## 15. Open questions

1. Test data: simulated DGPs for unit tests; which real dataset to
   cross-check against EViews?
2. ECM view: `adl_ecm(fit)` or a `form` argument.
3. Outlier critical values: fixed (3–3.5) or IIS-style levels tied to T.
4. Granger with several regressors: conditioning set and output format.

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
