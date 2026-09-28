#' Robust permutation test for equality of generalized variances
#'
#' Computes an MCD-based permutation test for
#' \eqn{H_0: |\Sigma_1| = \cdots = |\Sigma_k|} in a low-dimensional
#' multivariate setting. Group-specific robust shape transformations are
#' applied before the generalized-variance statistic is recalculated. The
#' scale component is retained; it is not divided out after transformation.
#'
#' @param x A numeric matrix or data frame. Rows are observations and columns
#'   are variables.
#' @param group A vector identifying the population or group of each row.
#' @param B Integer number of permutations. The default is 999.
#' @param alpha Numeric MCD coverage fraction in `(0.5, 1]`. The default is
#'   0.75.
#' @param seed Optional integer random seed.
#' @param na.rm Logical; if `TRUE`, rows containing missing values are removed.
#'   If `FALSE`, missing values produce an error.
#'
#' @return An object of class `RobPer_GVTest`, `MVTests`, and `list` with the
#'   following components: `statistic`, `p.value`, `successful`, `failed`,
#'   `permutation.statistics`, `B`, `alpha`, `group.sizes`, `p`, `k`, and
#'   `method`.
#'
#' @details
#' The permutation p-value is calculated as
#' \deqn{(1 + \sum_b I(T_b^* \geq T_{obs}))/(1+B_{ok}),}
#' where the sum is over successful permutations and `B_ok` is the number of
#' successful permutation fits. The observed and permuted statistics use the
#' same MCD-based transformation and refitting pipeline.
#'
#' If \eqn{S_i} is the group-specific MCD scatter matrix, the robust scale
#' component is \eqn{a_i = |S_i|^{1/p}} and the shape matrix is
#' \eqn{Q_i = S_i/a_i}. The transformed observations are multiplied by
#' \eqn{Q_i^{-1/2}}. No additional division by \eqn{\sqrt{a_i}} is applied,
#' because that would remove the generalized-variance information being tested.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("rrcov", quietly = TRUE)) {
#'   set.seed(1)
#'   x <- rbind(matrix(rnorm(60), 20, 3),
#'              matrix(rnorm(60), 20, 3))
#'   g <- rep(c("A", "B"), each = 20)
#'   fit <- RobPer_GVTest(x, g, B = 99, seed = 123)
#'   fit
#'   summary(fit)
#' }
#' }
#'
#' @export
RobPer_GVTest <- function(x, group, B = 999, alpha = 0.75,
                          seed = NULL, na.rm = TRUE) {
  if (!requireNamespace("rrcov", quietly = TRUE)) {
    stop("Package 'rrcov' is required for RobPer_GVTest.", call. = FALSE)
  }

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a numeric matrix or data frame.", call. = FALSE)
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("All columns of 'x' must be numeric.", call. = FALSE)
  }
  if (length(group) != nrow(x)) {
    stop("'group' must have one entry for each row of 'x'.", call. = FALSE)
  }
  if (anyNA(group)) {
    stop("'group' must not contain missing values.", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed)) {
      stop("'seed' must be a single finite numeric value.", call. = FALSE)
    }
    set.seed(seed)
  }
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) {
    stop("'B' must be a positive integer.", call. = FALSE)
  }
  B <- as.integer(B)
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0.5 || alpha > 1) {
    stop("'alpha' must be in the interval (0.5, 1].", call. = FALSE)
  }

  complete <- stats::complete.cases(x)
  if (any(!complete)) {
    if (!isTRUE(na.rm)) {
      stop("Missing values are present in 'x'.", call. = FALSE)
    }
    x <- x[complete, , drop = FALSE]
    group <- group[complete]
  }
  if (nrow(x) == 0L) {
    stop("No complete observations remain.", call. = FALSE)
  }

  group <- droplevels(factor(group))
  n_i <- tabulate(group)
  p <- ncol(x)
  k <- length(n_i)
  if (k < 2L) {
    stop("At least two groups are required.", call. = FALSE)
  }
  if (any(n_i <= p)) {
    stop("Each group must contain more observations than variables.",
         call. = FALSE)
  }

  logdet_spd <- function(S) {
    S <- (S + t(S)) / 2
    R <- tryCatch(chol(S), error = function(e) NULL)
    if (is.null(R)) return(NA_real_)
    2 * sum(log(diag(R)))
  }

  mcd_fit <- function(X) {
    fit <- rrcov::CovMcd(X, alpha = alpha)
    S <- (fit@cov + t(fit@cov)) / 2
    ld <- logdet_spd(S)
    if (!is.finite(ld)) {
      stop("MCD scatter is not positive definite.", call. = FALSE)
    }
    a <- exp(ld / ncol(X))
    list(center = as.numeric(fit@center), scatter = S,
         a = a, shape = S / a)
  }

  matrix_inv_sqrt <- function(S) {
    ee <- eigen((S + t(S)) / 2, symmetric = TRUE)
    if (any(!is.finite(ee$values)) || any(ee$values <= 0)) {
      stop("Shape matrix is not positive definite.", call. = FALSE)
    }
    ee$vectors %*% diag(1 / sqrt(ee$values), nrow = length(ee$values)) %*%
      t(ee$vectors)
  }

  gv_statistic <- function(scatter_list, sizes) {
    log_a <- vapply(scatter_list, function(S) {
      logdet_spd(S) / p
    }, numeric(1))
    if (any(!is.finite(log_a))) return(NA_real_)
    a <- exp(log_a)
    n <- sum(sizes)
    p * (n * log(sum(sizes * a) / n) - sum(sizes * log(a)))
  }

  x_list <- lapply(seq_len(k), function(i) {
    x[group == levels(group)[i], , drop = FALSE]
  })
  fits <- lapply(x_list, mcd_fit)

  # Remove group-specific shape/orientation while retaining robust scale.
  shape_adjusted <- lapply(seq_len(k), function(i) {
    fit <- fits[[i]]
    z <- sweep(x_list[[i]], 2, fit$center, "-")
    z %*% matrix_inv_sqrt(fit$shape)
  })

  # Refit transformed groups so the observed and permuted statistics share
  # exactly the same computational pipeline.
  adjusted_fits <- lapply(shape_adjusted, mcd_fit)
  T_obs <- gv_statistic(
    lapply(adjusted_fits, function(z) z$scatter), n_i
  )
  if (!is.finite(T_obs)) {
    stop("The observed robust statistic could not be computed.",
         call. = FALSE)
  }

  pooled <- do.call(rbind, shape_adjusted)
  T_perm <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    permuted <- pooled[sample.int(nrow(pooled)), , drop = FALSE]
    start <- 1L
    perm_groups <- lapply(n_i, function(size) {
      idx <- start:(start + size - 1L)
      start <<- start + size
      permuted[idx, , drop = FALSE]
    })
    perm_fits <- tryCatch(lapply(perm_groups, mcd_fit),
                          error = function(e) NULL)
    if (!is.null(perm_fits)) {
      T_perm[b] <- gv_statistic(
        lapply(perm_fits, function(z) z$scatter), n_i
      )
    }
  }

  ok <- is.finite(T_perm)
  if (!any(ok)) {
    stop("All MCD permutation fits failed.", call. = FALSE)
  }

  result <- list(
    statistic = T_obs,
    p.value = (1 + sum(T_perm[ok] >= T_obs)) / (1 + sum(ok)),
    method = "MCD permutation test for equality of generalized variances",
    data.name = deparse(substitute(x)),
    B = B,
    alpha = alpha,
    successful = sum(ok),
    failed = sum(!ok),
    permutation.statistics = T_perm,
    group.sizes = n_i,
    group.levels = levels(group),
    p = p,
    k = k,
    call = match.call()
  )
  class(result) <- c("RobPer_GVTest", "MVTests", "list")
  result
}

#' @export
print.RobPer_GVTest <- function(x, ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("Statistic:", format(x$statistic, digits = 6), "\n")
  cat("Permutation p-value:", format(x$p.value, digits = 6), "\n")
  cat("Groups:", x$k, " Variables:", x$p, "\n")
  cat("Successful permutations:", x$successful,
      " Failed:", x$failed, "\n")
  invisible(x)
}

#' @export
summary.RobPer_GVTest <- function(object, ...) {
  class(object) <- c("summary.RobPer_GVTest", "list")
  object
}

#' @export
print.summary.RobPer_GVTest <- function(x, ...) {
  print.RobPer_GVTest(x, ...)
  invisible(x)
}
