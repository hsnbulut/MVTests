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

The `RobCellT2_onesample()` function performs the cellMCD-based robust
one-sample Hotelling T2 test for testing a multivariate mean vector in the
presence of cellwise outliers. The constants required for the approximate
F distribution can be obtained using the `simRobCellT2_onesample()` function.

```r
library(MVTests)

set.seed(123)

X <- MASS::mvrnorm(
  n = 50,
  mu = rep(0, 5),
  Sigma = diag(5)
)

const <- simRobCellT2_onesample(
  n = 50,
  p = 5,
  nrep = 3000,
  seed = 123
)

fit <- RobCellT2_onesample(
  data = X,
  mu0 = rep(0, 5),
  d = const$d,
  q = const$q
)

summary(fit)
```
