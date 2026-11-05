# adltools (R)

[![R-CMD-check](https://github.com/maarranz/ADL-Tools-R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/maarranz/ADL-Tools-R/actions/workflows/R-CMD-check.yaml)

Teaching-oriented R tools for **autoregressive distributed lag (ADL)** models, developed for **Applied Time Series Econometrics (ATSE)** at **Universidad Carlos III de Madrid**.

> **Status: under development.** The package skeleton is in place; no functions are available yet. The design decisions are recorded in [`DESIGN.md`](DESIGN.md).

## Scope

`adltools` brings together, under one consistent interface, the tools needed to formulate, estimate and test ADL models as they are taught in the ATSE chapters on dynamic regression (ADL I–IV). Several of these tests are scattered across R packages, and some are not available elsewhere in the form used in the course.

The first version (0.1) will cover:

- **Formulation and estimation**: models written either with `dynlm`-style formulas (`y ~ L(y, 1) + L(x, 0:1)`) or with lag orders (`adl(y ~ x, data, order = c(1, 1))`), including the dead-start model and unlagged regressors such as outlier dummies.
- **Inference**: coefficient tables in the chapters' α/β notation, joint Wald tests and the long-run multiplier with delta-method standard errors, with conventional or Newey–West (HAC) covariance.
- **Granger non-causality**: the direct (Sargent) procedure of the chapters, with separate tests for simple causality in each direction and for instantaneous causality, and the classification of the relation into one of the eight causal types.

Later versions will add diagnostic and stability tests, outlier detection, dynamic effects and dynamic forecasting (see `DESIGN.md`).

## Related packages

`adltools` is part of a suite of teaching packages for ATSE:

| Package | Repository | Covers |
|:--|:--|:--|
| arimatools | [ARIMA-Tools-R](https://github.com/maarranz/ARIMA-Tools-R) | ARIMA estimation, diagnostics, outliers and dynamic forecasting |
| forecastingtools | [Forecasting-Tools-R](https://github.com/maarranz/Forecasting-Tools-R) | Static forecasting, forecast accuracy, evaluation and comparison (including ADL forecasts) |
| unitroottools | [UnitRoot-Tools-R](https://github.com/maarranz/UnitRoot-Tools-R) | Unit-root and stationarity testing |
| adltools | ADL-Tools-R | ADL models (this package) |

## Installation

Once the first version is released:

```r
# install.packages("remotes")
remotes::install_github("maarranz/ADL-Tools-R")
library(adltools)
```

Requires R ≥ 4.1.

## Author

**Miguel A. Arranz**
Universidad Carlos III de Madrid
<maarranz@eco.uc3m.es>
ORCID: [0000-0002-5951-2284](https://orcid.org/0000-0002-5951-2284)

## Licence

MIT © Miguel A. Arranz
