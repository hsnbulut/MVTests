# **MVTests: Multivariate Hypothesis Tests**

- I have constructed this package to be able to use parametric and non-parametric multivariate hypothesis tests in R. 
- I will develop this package continuously.

# To install the package from GitHub

## Install devtools if needed
install.packages("devtools")

## Install MVTests from GitHub
devtools::install_github("hsnbulut/MVTests")
library(MVTests)

# Examples

## Adaptive Wrapped Robust Canonical Correlation Analysis

`AWRcca()` combines cellwise wrapping, analytical shrinkage, and MCD
reweighting in an initial canonical-score space. The binary weights enter a
second shrinkage-regularized CCA fit, so the estimator can be used when the
combined dimension is larger than the sample size.

```r
library(MVTests)

set.seed(123)
n <- 60
p <- 40
q <- 50
latent <- rnorm(n)
X <- outer(latent, seq(0.80, 0.30, length.out = p)) +
  matrix(rnorm(n * p), n, p)
Y <- outer(latent, seq(0.75, 0.25, length.out = q)) +
  matrix(rnorm(n * q), n, q)

fit <- AWRcca(X, Y)
fit$cor[1]
fit$retained
fit$shrink_final
```



## Robust High-Dimensional MANOVA Test

The `RobHDMANOVA()` function performs a weighted MRCD-based robust MANOVA test for high-dimensional data.

```r
library(MVTests)

fit <- RobHDMANOVA(x = X, group = group, N = 999,
                   alpha = 0.75, tau = 0.975,
                   cutoff = "normal")

summary(fit)
```

## Cellwise Robust Two-Sample Hotelling T2 Test

The `CellMCDT2()` function performs the cellMCD-based robust permutation Hotelling T2 test for comparing the mean vectors of two independent groups in the presence of cellwise outliers.

```r
library(MVTests)

set.seed(123)
X1 <- mvtnorm::rmvnorm(n = 30, mean = rep(0, 5), sigma = diag(5))
X2 <- mvtnorm::rmvnorm(n = 30, mean = rep(0, 5), sigma = diag(5))

fit <- CellMCDT2(X1 = X1, X2 = X2, B = 199, seed = 123)
summary(fit)
```

## Cellwise Robust One-Sample Hotelling T2 Test

`RobCellT2_onesample()` tests a multivariate mean with the cellMCD location and
scatter estimates. In version 2.3.6, calibrate its empirical scaled-F reference
with the following interface:

```r
simRobCellT2_onesample(X = NULL, mu0, Sigma = NULL, n = NULL, p = NULL, nper = 3000)
```

`mu0` is required and remains the simulation mean. With X, n and p come from its
rows and columns; explicitly supplied dimensions must agree. If Sigma is
supplied, it is used directly. Otherwise X is required, and the bias-corrected
scatter S from `cellWise::cellMCD(X)` on **raw X** supplies the covariance. The
observed location estimate never replaces mu0. Without X, n, p and Sigma are all
required. X must be finite and complete: remove incomplete rows before calibration
and use the same observations in the test.

The covariance fit and every simulated fit use default cellMCD settings
(alpha=.75, quant=.99, crit=1e-4). Every simulation sample comes from
N_p(mu0,Sigma), and T2 is measured against the fixed mu0. The returned d,q match
Monte Carlo moments and require finite q>4. Returned metadata includes moments,
Sigma, mu0, n, p, nper, successful/failed counts, failure messages and fit settings.
Use `set.seed()` for reproducibility. A small nper can yield an unstable moment
match; the function stops instead of returning unusable constants.

This is an empirical covariance-dependent approximation, not an exact F law.
A covariance estimated from X is a plug-in estimate, not guaranteed to equal
population covariance. Standardization does not remove dependence on correlations.
General whitening changes the coordinate-wise cell structure and is not routine
preprocessing. Match n, p, covariance and cellMCD settings when applying the test.

```r
library(MVTests)
mu0 <- c(1, -1)
Sigma_true <- matrix(c(2, 0.4, 0.4, 1), 2)
set.seed(123)
X <- MASS::mvrnorm(n = 30, mu = mu0, Sigma = Sigma_true)

set.seed(1)
const <- simRobCellT2_onesample(X = X, mu0 = mu0, nper = 50)
fit <- RobCellT2_onesample(data = X, mu0 = mu0, d = const$d, q = const$q)
fit$p.value

# Calibration under an explicitly specified covariance:
set.seed(1)
const_explicit <- simRobCellT2_onesample(
  mu0 = mu0, Sigma = Sigma_true, n = 30, p = 2, nper = 50
)
```

For production calibration use nper=3000; the smaller example count limits runtime.

### Migration from the previous calibration API

The six-argument interface does not accept nrep, seed, data or additional cellMCD
arguments. Use nper for nrep, X for data, explicit mu0, and set.seed() for seed.
Old positional calls must not be interpreted as new X/mu0 calls. The exported
`simRobCellT2_onesample_legacy()` function preserves the old positional argument
order, custom cellMCD settings, default identity covariance warning, and the
previous centered/fixed-center data-based calibration design. It emits a
deprecation warning. This provides an explicit migration path; the new raw-X
covariance estimate differs statistically from the legacy fixed-center estimate.

```r
# Existing identity-model call, using the compatibility interface:
legacy <- suppressWarnings(simRobCellT2_onesample_legacy(
  n = 30, p = 2, nrep = 30, seed = 1, Sigma = diag(2)
))
# Corresponding new call (explicit null mean/covariance):
set.seed(1)
current <- simRobCellT2_onesample(mu0 = rep(0, 2), Sigma = diag(2),
                                n = 30, p = 2, nper = 30)
```
