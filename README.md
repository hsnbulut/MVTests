# **MVTests: Multivariate Hypothesis Tests**

- I have constructed this package to be able to use parametric and non-parametric multivariate hypothesis tests in R. 
- I will develop this package continuously.

# To install the package from GitHub

## Install devtools if needed
install.packages("devtools")

## Install MVTests from GitHub
devtools::install_github("hsnbulut/MVTests")
library(MVTests)


## Robust High-Dimensional MANOVA Test

The `RobHDMANOVA()` function performs a weighted MRCD-based robust MANOVA test for high-dimensional data.

```r
library(MVTests)

fit <- RobHDMANOVA(x = X, group = group, N = 999,
                   alpha = 0.75, tau = 0.975,
                   cutoff = "normal")

summary(fit)

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
