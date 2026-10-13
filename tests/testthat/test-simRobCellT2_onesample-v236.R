skip_if_not_installed("MASS")
skip_if_not_installed("cellWise")

test_that("2.3.6 signature is exact and mu0 is required", {
  expect_equal(names(formals(simRobCellT2_onesample)),
               c("X", "mu0", "Sigma", "n", "p", "nper"))
  expect_error(simRobCellT2_onesample(), "mu0 is required")
  expect_error(simRobCellT2_onesample(X = matrix(1, 30, 2)), "mu0 is required")
  expect_error(simRobCellT2_onesample(30, 2), "legacy")
})

test_that("raw X supplies its covariance and dimensions with reproducible results", {
  mu0 <- c(1, -1)
  set.seed(123)
  X <- MASS::mvrnorm(30, mu0, matrix(c(2, .4, .4, 1), 2))
  colnames(X) <- c("x", "y")
  set.seed(1)
  result <- simRobCellT2_onesample(as.data.frame(X), mu0, nper = 30)
  expect_equal(c(result$n, result$p, result$nper), c(30, 2, 30))
  expect_equal(result$mu0, mu0)
  expect_equal(result$mode, "data")
  expect_equal(result$covariance.source, "cellMCD S from raw X")
  expect_equal(result$n.success + result$n.failed, 30)
  expect_equal(result$cellMCD.settings, result$covariance.settings)
  expect_false(result$covariance.settings$fixedCenter)
  reference <- NULL
  set.seed(1)
  invisible(capture.output(reference <- cellWise::cellMCD(X)))
  expect_equal(result$Sigma, reference$S)
  set.seed(1)
  expect_equal(result, simRobCellT2_onesample(X, mu0, n = 30, p = 2, nper = 30))
})

test_that("explicit Sigma overrides X estimation and matches the no-X mode", {
  Sigma <- matrix(c(2, .4, .4, 1), 2)
  X <- matrix(0, 30, 2) # cellMCD would fail on this constant matrix
  set.seed(1)
  with.X <- simRobCellT2_onesample(X, c(1, -1), Sigma = Sigma, nper = 30)
  set.seed(1)
  without.X <- simRobCellT2_onesample(mu0 = c(1, -1), Sigma = Sigma,
                                      n = 30, p = 2, nper = 30)
  expect_equal(with.X, without.X)
  expect_equal(without.X$Sigma, Sigma)
  expect_equal(without.X$mode, "explicit")
  expect_null(without.X$covariance.settings)
})

test_that("simulation uses fixed mu0, supplied covariance and test T2 calculation", {
  mu0 <- c(7, -3)
  Sigma <- matrix(c(2, .4, .4, 1), 2)
  raw <- matrix(seq_len(60), 30, 2)
  values <- c(rep(1, 9), 25)
  calls <- list(); generated <- list()
  local_mocked_bindings(cellMCD = function(X, fixedCenter = FALSE, ...) {
    calls[[length(calls) + 1L]] <<- X
    if (length(calls) == 1L) return(list(mu = c(100, 200), S = Sigma))
    k <- length(calls) - 1L
    list(mu = mu0 + c(sqrt(values[k] / nrow(X)), 0), S = diag(2))
  }, .package = "cellWise")
  local_mocked_bindings(mvrnorm = function(n, mu, Sigma, ...) {
    generated[[length(generated) + 1L]] <<- list(n = n, mu = mu, Sigma = Sigma)
    matrix(rep(mu, each = n), n, length(mu))
  }, .package = "MASS")
  result <- simRobCellT2_onesample(raw, mu0, nper = 10)
  expect_equal(calls[[1]], raw) # no centering or fitted-location substitution
  expect_equal(result$Sigma, Sigma)
  expect_equal(length(generated), 10)
  for (sample in generated) expect_equal(sample, list(n = 30L, mu = mu0, Sigma = Sigma))
  expect_equal(result$mean.T2, mean(values))
  expect_equal(result$var.T2, var(values))
  q <- 4 / (var(values) / mean(values)^2 - 1) + 4
  expect_equal(result$q, q)
  expect_equal(result$d, mean(values) * (q - 2) / q)
})

test_that("2.3.6 validates dimensions, null means and covariance", {
  X <- matrix(1, 30, 2)
  expect_error(simRobCellT2_onesample(X, c(0, 0), n = 29), "equal nrow")
  expect_error(simRobCellT2_onesample(X, c(0, 0), p = 3), "equal ncol")
  for (mu in list(1, c(NA, 0), c(Inf, 0), c("a", "b"), matrix(1, 1, 2))) {
    expect_error(simRobCellT2_onesample(X, mu), "finite numeric vector of length p")
  }
  expect_error(simRobCellT2_onesample(matrix("a", 30, 2), c(0, 0)), "finite numeric")
  expect_error(simRobCellT2_onesample(matrix(NA_real_, 30, 2), c(0, 0)), "missing values")
  expect_error(simRobCellT2_onesample(mu0 = c(0, 0), n = 30, p = 2), "Without X")
  expect_error(simRobCellT2_onesample(mu0 = c(0, 0), Sigma = diag(2)), "Without X")
  for (number in list(NA_real_, Inf, 2.5, 1)) {
    expect_error(simRobCellT2_onesample(mu0 = c(0, 0), n = number, p = 2,
                                       Sigma = diag(2)), "n must be a finite integer")
    expect_error(simRobCellT2_onesample(mu0 = c(0, 0), n = 30, p = number,
                                       Sigma = diag(2)), "p must be a finite integer")
  }
  for (number in list(NA_real_, Inf, .5, 0)) {
    expect_error(simRobCellT2_onesample(X, c(0, 0), Sigma = diag(2), nper = number),
                 "nper must be a positive finite integer")
  }
  for (S in list(diag(3), matrix(NA_real_, 2, 2), matrix("a", 2, 2))) {
    expect_error(simRobCellT2_onesample(X, c(0, 0), Sigma = S), "finite numeric p-by-p")
  }
  expect_error(simRobCellT2_onesample(X, c(0, 0), Sigma = matrix(c(1, 0, 1, 1), 2)),
               "symmetric")
  expect_error(simRobCellT2_onesample(X, c(0, 0), Sigma = matrix(1, 2, 2)), "positive definite")
  expect_error(simRobCellT2_onesample(X, c(0, 0), Sigma = diag(c(1, -1))), "positive definite")
})

test_that("raw-X covariance failures and invalid scatter produce clear errors", {
  X <- matrix(seq_len(60), 30, 2)
  local_mocked_bindings(cellMCD = function(X, ...) stop("failed fit"), .package = "cellWise")
  expect_error(simRobCellT2_onesample(X, c(0, 0)), "raw X failed: failed fit")
  local_mocked_bindings(cellMCD = function(X, ...) list(raw.S = diag(2)), .package = "cellWise")
  expect_error(simRobCellT2_onesample(X, c(0, 0)), "bias-corrected scatter S")
  local_mocked_bindings(cellMCD = function(X, ...) list(S = diag(c(1, -1))), .package = "cellWise")
  expect_error(simRobCellT2_onesample(X, c(0, 0)), "cellMCD covariance S must be positive definite")
})

test_that("replication failures are counted and retain informative diagnostics", {
  i <- 0L
  values <- c(rep(1, 9), 25)
  local_mocked_bindings(cellMCD = function(X, ...) {
    i <<- i + 1L
    if (i <= 2L) stop("replicate fit failed")
    list(mu = c(sqrt(values[i - 2L] / nrow(X)), 0), S = diag(2))
  }, .package = "cellWise")
  result <- simRobCellT2_onesample(mu0 = c(0, 0), Sigma = diag(2), n = 30, p = 2, nper = 12)
  expect_equal(result$n.success, 10)
  expect_equal(result$n.failed, 2)
  expect_equal(names(result$failure.messages), c("1", "2"))
  expect_true(all(grepl("replicate fit failed", result$failure.messages)))
  local_mocked_bindings(cellMCD = function(X, ...) stop("all failed"), .package = "cellWise")
  expect_error(simRobCellT2_onesample(mu0 = c(0, 0), Sigma = diag(2), n = 30, p = 2, nper = 10),
               "Too few successful.*0/10")
})

test_that("unstable F moment matches fail explicitly", {
  local_mocked_bindings(cellMCD = function(X, ...) list(mu = c(1, 1), S = diag(2)),
                        .package = "cellWise")
  expect_error(simRobCellT2_onesample(mu0 = c(0, 0), Sigma = diag(2), n = 30, p = 2, nper = 10),
               "moments must be positive and finite")
  i <- 0L
  local_mocked_bindings(cellMCD = function(X, ...) {
    i <<- i + 1L
    list(mu = c(sqrt(i / nrow(X)), 0), S = diag(2))
  }, .package = "cellWise")
  expect_error(simRobCellT2_onesample(mu0 = c(0, 0), Sigma = diag(2), n = 30, p = 2, nper = 10),
               "finite-variance F distribution.*q > 4")
})

test_that("explicit legacy migration preserves old positional results", {
  expect_warning(old <- simRobCellT2_onesample_legacy(30, 2, 30, .75, .99, 1e-4, 1, diag(2)),
                 "legacy calibration API is deprecated")
  set.seed(1)
  new <- simRobCellT2_onesample(mu0 = c(0, 0), Sigma = diag(2), n = 30, p = 2, nper = 30)
  expect_equal(old[c("d", "q", "mean.T2", "var.T2", "Sigma")],
               new[c("d", "q", "mean.T2", "var.T2", "Sigma")])
})
