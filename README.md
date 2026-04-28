# **MVTests: Multivariate Hypothesis Tests**

- I have constructed this package to be able to use parametric and non-parametric multivariate hypothesis tests in R. 
- I will develop this package continuously.

# To install the package from GitHub

## Install devtools if needed
install.packages("devtools")

## Install MVTests from GitHub
devtools::install_github("hsnbulut/MVTests")
library(MVTests)


## Weighted MRCD-Based Robust MANOVA Test

The `RobHDMANOVA()` function performs a weighted MRCD-based robust MANOVA test for high-dimensional data.

```r
library(MVTests)

fit <- RobHDMANOVA(x = X, group = group, N = 999,
                   alpha = 0.75, tau = 0.975,
                   cutoff = "normal")

summary(fit)
