skip_if_not_installed("MASS")
skip_if_not_installed("cellWise")

test_that("data-based calibration estimates S and is reproducible", {
  set.seed(123)
  mu0 <- c(1, -1)
  X <- MASS::mvrnorm(30, mu0, matrix(c(2, 0.4, 0.4, 1), 2))
  colnames(X) <- c("x", "y")
  result <- .simRobCellT2_onesample_legacy_impl(data = as.data.frame(X), mu0 = mu0,
                                  nrep = 30, seed = 1)
  expect_equal(result$mode, "data")
  expect_equal(c(result$n, result$p), c(30, 2))
  expect_equal(result$mu0, mu0)
  expect_equal(result$covariance.source, "cellMCD S")
  expect_true(result$covariance.settings$fixedCenter)
  expect_false(result$cellMCD.settings$fixedCenter)
  reference <- NULL
  set.seed(1)
  invisible(capture.output({
    reference <- do.call(cellWise::cellMCD,
      c(list(X = sweep(X, 2, mu0, "-")), result$covariance.settings))
  }))
  expect_equal(result$Sigma, reference$S)
  expect_equal(result, .simRobCellT2_onesample_legacy_impl(data = X, mu0 = mu0,
                                             nrep = 30, seed = 1))
  expect_equal(result$n.success + result$n.failed, result$nrep)
})

test_that("complete-case handling agrees with the test", {
  set.seed(123)
  X <- MASS::mvrnorm(30, c(1, -1), matrix(c(2, 0.4, 0.4, 1), 2))
  incomplete <- rbind(X, c(NA, 0), c(1, NaN))
  result <- .simRobCellT2_onesample_legacy_impl(data = incomplete, mu0 = c(1, -1),
                                  nrep = 30, seed = 1)
  complete <- .simRobCellT2_onesample_legacy_impl(data = X, mu0 = c(1, -1), nrep = 30,
                                    seed = 1, na.rm = FALSE)
  expect_equal(result$n, 30)
  expect_equal(result$n.omitted, 2)
  expect_equal(result[c("d", "q", "Sigma")], complete[c("d", "q", "Sigma")])
  test <- RobCellT2_onesample(incomplete, c(1, -1), result$d, result$q)
  expect_equal(unname(test$parameter[c("n", "p")]), c(result$n, result$p))
  missing.cells <- X
  missing.cells[1, 1] <- NA_real_
  retained <- .simRobCellT2_onesample_legacy_impl(data = missing.cells, mu0 = c(1, -1),
                                    nrep = 30, seed = 1, na.rm = FALSE)
  retained.test <- RobCellT2_onesample(missing.cells, c(1, -1),
                                      retained$d, retained$q, na.rm = FALSE)
  expect_equal(retained$n.omitted, 0)
  expect_equal(unname(retained.test$parameter[c("n", "p")]),
               c(retained$n, retained$p))
  reference <- NULL
  set.seed(1)
  invisible(capture.output({
    reference <- do.call(cellWise::cellMCD,
      c(list(X = sweep(missing.cells, 2, c(1, -1), "-")),
        retained$covariance.settings))
  }))
  expect_equal(retained$Sigma, reference$S)
})

test_that("samples use fixed mu0 and plug-in S, with shared settings", {
  mu0 <- c(7, -3)
  X <- matrix(seq_len(60) / 10, 30, 2)
  Sigma <- matrix(c(2, 0.4, 0.4, 1), 2)
  calls <- list()
  generated <- list()
  values <- c(rep(1, 9), 25)
  local_mocked_bindings(cellMCD = function(X, alpha = 0.75, quant = 0.99,
      crit = 1e-04, noCits = 100, fixedCenter = FALSE, ...) {
    calls[[length(calls) + 1L]] <<- list(X = X, alpha = alpha, quant = quant,
      crit = crit, noCits = noCits, fixedCenter = fixedCenter)
    if (fixedCenter) return(list(mu = c(100, 100), S = Sigma, raw.S = diag(2)))
    k <- length(calls) - 1L
    list(mu = mu0 + c(sqrt(values[k] / nrow(X)), 0), S = diag(2))
  }, .package = "cellWise")
  local_mocked_bindings(mvrnorm = function(n, mu, Sigma, ...) {
    generated[[length(generated) + 1L]] <<- list(n = n, mu = mu, Sigma = Sigma)
    matrix(rep(mu, each = n), n, length(mu))
  }, .package = "MASS")
  result <- .simRobCellT2_onesample_legacy_impl(data = X, mu0 = mu0, nrep = 10, seed = 1,
                                  alpha = 0.8, quant = 0.98, crit = 1e-05,
                                  noCits = 20, fixedCenter = FALSE)
  expect_equal(calls[[1]]$X, sweep(X, 2, mu0, "-"))
  expect_equal(result$Sigma, Sigma)
  expect_equal(length(generated), 10)
  for (sample in generated) {
    expect_equal(sample, list(n = 30L, mu = mu0, Sigma = Sigma))
  }
  for (call in calls) {
    expect_equal(call[c("alpha", "quant", "crit", "noCits")],
                 list(alpha = 0.8, quant = 0.98, crit = 1e-05, noCits = 20))
  }
  expect_true(calls[[1]]$fixedCenter)
  expect_true(all(vapply(calls[-1], function(x) !x$fixedCenter, logical(1))))
  expect_equal(result$mean.T2, mean(values))
  expect_equal(result$var.T2, var(values))
  q <- (2 + 2) / (var(values) / mean(values)^2 - 1) + 4
  expect_equal(result$q, q)
  expect_equal(result$d, mean(values) * (q - 2) / q)
})

test_that("malformed inputs and conflicting modes fail clearly", {
  X <- matrix(seq_len(60), 30, 2)
  for (mu in list(1, c(NA, 1), c(Inf, 1), c("a", "b"), matrix(1, 1, 2))) {
    expect_error(.simRobCellT2_onesample_legacy_impl(data = X, mu0 = mu),
                 "finite numeric vector of length p")
  }
  expect_error(.simRobCellT2_onesample_legacy_impl(data = X), "mu0 is required")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = matrix("a", 30, 2), mu0 = c(0, 0)),
               "numeric matrix or data frame")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = matrix(Inf, 30, 2), mu0 = c(0, 0)),
               "finite values")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = matrix(NA_real_, 30, 2),
                                     mu0 = c(0, 0)), "n must be a finite integer")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = matrix(1, 30, 1), mu0 = 0),
               "p must be a finite integer")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = X, mu0 = c(0, 0), na.rm = NA),
               "na.rm must be TRUE or FALSE")
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, data = X, mu0 = c(0, 0)),
               "not both modes")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = X, mu0 = c(0, 0), Sigma = diag(2)),
               "not both modes")
  for (value in list(NA_real_, Inf, 1.5, 0)) {
    expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, nrep = value, Sigma = diag(2)),
                 "nrep must be a positive finite integer")
  }
  for (value in list(NA_real_, Inf, 2.5, 1)) {
    expect_error(.simRobCellT2_onesample_legacy_impl(value, 2, Sigma = diag(2)),
                 "n must be a finite integer")
    expect_error(.simRobCellT2_onesample_legacy_impl(30, value, Sigma = diag(2)),
                 "p must be a finite integer")
  }
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, mu0 = c(0, Inf), Sigma = diag(2)),
               "finite numeric vector")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = X, mu0 = c(0, 0),
                                     fixedCenter = TRUE), "free replicate locations")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = X, mu0 = c(0, 0),
                                     fixedCenter = 1), "free replicate locations")
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, Sigma = diag(2), unknownSetting = TRUE),
               "Invalid cellMCD arguments")
  expect_error(.simRobCellT2_onesample_legacy_impl(30, 2, Sigma = diag(2), X = X),
               "Invalid cellMCD arguments")
})

test_that("failed plug-in fits never fall back to classical covariance", {
  X <- matrix(seq_len(60), 30, 2)
  local_mocked_bindings(cellMCD = function(X, fixedCenter = FALSE, ...) {
    stop("fit failed")
  }, .package = "cellWise")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = X, mu0 = c(0, 0)),
               "plug-in null covariance estimation failed: fit failed")
})

test_that("invalid plug-in scatter is rejected", {
  X <- matrix(seq_len(60), 30, 2)
  for (S in list(NULL, diag(3), matrix(NA_real_, 2, 2),
                 matrix(c(1, 0, 1, 1), 2), diag(c(1, -1)))) {
    local_mocked_bindings(cellMCD = function(X, fixedCenter = FALSE, ...) {
      list(S = S)
    }, .package = "cellWise")
    expect_error(.simRobCellT2_onesample_legacy_impl(data = X, mu0 = c(0, 0)),
                 "cellMCD plug-in covariance S must be")
  }
})

test_that("insufficient or invalid replicate fits fail clearly", {
  X <- matrix(seq_len(60), 30, 2)
  for (replicate.fit in list(NULL, list(mu = c(Inf, 0), S = diag(2)),
                             list(mu = c(0, 0), S = diag(2)),
                             list(mu = c(1, 0), S = diag(c(1, -1))))) {
    local_mocked_bindings(cellMCD = function(X, fixedCenter = FALSE, ...) {
      if (fixedCenter) return(list(S = diag(2)))
      if (is.null(replicate.fit)) stop("replicate failed")
      replicate.fit
    }, .package = "cellWise")
    expect_error(.simRobCellT2_onesample_legacy_impl(data = X, mu0 = c(0, 0), nrep = 10),
                 "Too few successful cellMCD replications")
  }
})

test_that("data-based mode also rejects unstable moments", {
  local_mocked_bindings(cellMCD = function(X, fixedCenter = FALSE, ...) {
    if (fixedCenter) return(list(S = diag(2)))
    list(mu = c(1, 1), S = diag(2))
  }, .package = "cellWise")
  expect_error(.simRobCellT2_onesample_legacy_impl(data = matrix(seq_len(60), 30, 2),
                                     mu0 = c(0, 0), nrep = 10),
               "moments must be positive and finite")
})

test_that("older APIs warn and record the free-center plug-in exception", {
  i <- 0L
  values <- c(rep(1, 9), 25)
  local_mocked_bindings(cellMCD = function(X, ...) {
    i <<- i + 1L
    if (i == 1L) return(list(S = diag(2)))
    list(mu = c(sqrt(values[i - 1L] / nrow(X)), 0), S = diag(2))
  }, .package = "cellWise")
  expect_warning(
    result <- .simRobCellT2_onesample_legacy_impl(data = matrix(seq_len(60), 30, 2),
                                    mu0 = c(0, 0), nrep = 10, seed = 1),
    "does not support fixedCenter"
  )
  expect_false("fixedCenter" %in% names(result$covariance.settings))
  expect_equal(result$covariance.settings, result$cellMCD.settings)
})
