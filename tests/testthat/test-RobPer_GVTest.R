test_that("RobPer_GVTest returns the expected object", {
  skip_if_not_installed("rrcov")

  data(crabs, package = "MASS")
  vars <- c("FL", "RW", "CL", "CW", "BD")
  keep <- paste0(crabs$sp, crabs$sex) %in% c("BM", "BF", "OM")
  x <- as.matrix(crabs[keep, vars])
  group <- paste0(crabs$sp[keep], crabs$sex[keep])

  fit <- RobPer_GVTest(
    x = x,
    group = group,
    B = 19,
    alpha = 0.75,
    seed = 20260927
  )

  expect_s3_class(fit, "RobPer_GVTest")
  expect_s3_class(fit, "MVTests")
  expect_true(is.numeric(fit$statistic))
  expect_true(is.numeric(fit$p.value))
  expect_length(fit$permutation.statistics, 19)
  expect_equal(fit$successful + fit$failed, 19)
  expect_true(fit$successful > 0)
})

test_that("RobPer_GVTest is reproducible for a fixed seed", {
  skip_if_not_installed("rrcov")

  set.seed(11)
  x <- rbind(
    matrix(rnorm(60), nrow = 20, ncol = 3),
    matrix(rnorm(60), nrow = 20, ncol = 3)
  )
  group <- rep(c("A", "B"), each = 20)

  fit1 <- RobPer_GVTest(x, group, B = 19, seed = 101)
  fit2 <- RobPer_GVTest(x, group, B = 19, seed = 101)

  expect_equal(fit1$statistic, fit2$statistic)
  expect_equal(fit1$p.value, fit2$p.value)
  expect_equal(fit1$permutation.statistics,
               fit2$permutation.statistics)
})

test_that("RobPer_GVTest validates the number of groups", {
  skip_if_not_installed("rrcov")

  x <- matrix(rnorm(60), nrow = 20, ncol = 3)
  expect_error(
    RobPer_GVTest(x, rep("A", 20), B = 9),
    "At least two groups are required"
  )
})
