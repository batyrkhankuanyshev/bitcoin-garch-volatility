# Bitcoin Volatility Analysis — GARCH(1,1)

**EF4822 Group 1 (Thursday)**

Analysis of Bitcoin daily return volatility using a GARCH(1,1) model with Student-t innovations, covering the period **2021–2024**.

## Overview

| Item | Detail |
|------|--------|
| Asset | Bitcoin (BTC/USD) |
| Data source | CoinMarketCap via `crypto2` R package |
| Period | 2021-01-01 → 2024-12-31 |
| Model | GARCH(1,1) with Student-t distribution |
| Benchmark | Constant variance model |

## Files

| File | Description |
|------|-------------|
| `bitcoin_garch_analysis.R` | Full R script: data download, EDA, GARCH estimation, diagnostics |
| `EF4822_Project_Brief.pdf` | Project assignment brief |

## What the Script Does

1. **Data** — Downloads daily BTC prices from CoinMarketCap and computes log returns
2. **Descriptive stats** — Mean, std dev, skewness, excess kurtosis
3. **Plots** — Price series, returns, squared returns, ACF/PACF
4. **ARCH test** — Ljung-Box test on squared returns
5. **GARCH(1,1)** — Estimates ω, α, β with robust standard errors
6. **Model comparison** — GARCH vs constant variance (AIC, BIC, log-likelihood)
7. **Diagnostics** — Residual tests, QQ-plot, news impact curve
8. **Output** — Saves results to CSV and TXT files

## Requirements

```r
install.packages(c("crypto2", "moments", "rugarch", "tseries"))
```

## Usage

```r
source("bitcoin_garch_analysis.R")
```
