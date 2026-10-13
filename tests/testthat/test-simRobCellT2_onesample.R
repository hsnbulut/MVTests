skip_if_not_installed("MASS")
skip_if_not_installed("cellWise")

# Independent reference to the original identity-model simulation and moments.
legacy_identity_calibration <- function(n, p, nrep, seed) {
  set.seed(seed)
  values <- replicate(nrep, {
    X <- MASS::mvrnorm(n, rep(0, p), diag(p))
    fit <- NULL
    invisible(capture.output({
      fit <- tryCatch(cellWise::cellMCD(X, alpha = 0.75, quant = 0.99,
                                      crit = 1e-04), error = function(e) NULL)
    }))
    if (is.null(fit)) return(NA_real_)
    mu <- if (!is.null(fit$mu)) fit$mu else fit$center
    S <- if (!is.null(fit$S)) fit$S else fit$cov
    if (is.null(mu) || is.null(S)) return(NA_real_)
    delta <- as.numeric(mu)
    tryCatch(as.numeric(n * t(delta) %*% solve(as.matrix(S)) %*% delta),
             error = function(e) NA_real_)
  })
  values <- values[is.finite(values) & values > 0]
  m <- mean(values)
  v <- var(values)
  q <- (p + 2) / ((v / m^2) * (p / 2) - 1) + 4
  list(d = m * (q - 2) / q, q = q, mean.T2 = m, var.T2 = v,
       n.success = length(values))
}

test_that("explicit identity Sigma reproduces legacy calibration", {
  legacy <- legacy_identity_calibration(30, 2, 30, 1)
  result <- .simRobCellT2_onesample_legacy_impl(30, 2, 30, 0.75, 0.99, 1e-04, 1,
                                  Sigma = diag(2))
  expect_equal(result[names(legacy)], legacy)
  expect_equal(names(formals(.simRobCellT2_onesample_legacy_impl)),
               c("n", "p", "nrep", "alpha", "quant", "crit", "seed",
                 "Sigma", "...", "data", "mu0", "na.rm"))
})

test_that("correlated positive-definite Sigma is accepted and returned", {
  Sigma <- matrix(c(2, 0.4, 0.4, 1), 2)
  result <- .simRobCellT2_onesample_legacy_impl(30, 2, 30, seed = 1, Sigma = Sigma)
  expect_equal(result$Sigma, Sigma)
  expect_true(is.finite(result$d) && result$d > 0)
  expect_true(is.finite(result$q) && result$q > 4)
})

test_that("Sigma validation gives informative errors", {
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, Sigma = diag(3)),
               "finite numeric p-by-p matrix")
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, Sigma = matrix(NA_real_, 2, 2)),
               "finite numeric p-by-p matrix")
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, Sigma = matrix("a", 2, 2)),
               "finite numeric p-by-p matrix")
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, Sigma = matrix(c(1, 0, 1, 1), 2)),
               "symmetric")
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, Sigma = matrix(1, 2, 2)),
               "positive definite")
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, Sigma = diag(c(1, -1))),
               "positive definite")
})

test_that("omitted Sigma warns and old positional calls retain their meaning", {
  expect_warning(
    result <- .simRobCellT2_onesample_legacy_impl(30, 2, 30, 0.75, 0.99, 1e-04, 1),
    "independence/identity-covariance model.*correlated data"
  )
  explicit <- .simRobCellT2_onesample_legacy_impl(30, 2, 30, 0.75, 0.99, 1e-04, 1, diag(2))
  expect_equal(result$covariance.source, "default identity")
  expect_equal(explicit$covariance.source, "supplied Sigma")
  shared <- setdiff(names(result), "covariance.source")
  expect_equal(result[shared], explicit[shared])
})

test_that("incompatible moments do not return unusable constants", {
  # Deterministic fit isolates the finite-variance guard from cellMCD randomness.
  local_mocked_bindings(cellMCD = function(X, ...) {
    list(mu = rep(1, ncol(X)), S = diag(ncol(X)))
  }, .package = "cellWise")
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, 10, seed = 1, Sigma = diag(2)),
               "moments must be positive and finite")
})

test_that("positive but incompatible or near-boundary moments are rejected", {
  for (values in list(seq_len(10),
                      c(rep(1, 9), (9 + sqrt(10)) / (sqrt(10) - 1)))) {
    i <- 0L
    local_mocked_bindings(cellMCD = function(X, ...) {
      i <<- i + 1L
      list(mu = c(sqrt(values[i] / nrow(X)), 0), S = diag(2))
    }, .package = "cellWise")
    expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, 10, seed = 1, Sigma = diag(2)),
                 "finite-variance F distribution.*q > 4")
  }
})
