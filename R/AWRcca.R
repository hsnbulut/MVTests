#' Adaptive Wrapped Robust Canonical Correlation Analysis
#'
#' Fits the Adaptive Wrapped Robust Canonical Correlation Analysis (AWRcca)
#' estimator. The procedure combines columnwise robust wrapping, analytical
#' shrinkage of the joint correlation matrix, and MCD reweighting in a
#' low-dimensional canonical score space. The score-space weights enter a
#' second shrinkage-regularized CCA fit, which remains well defined when
#' \eqn{p+q>n}.
#'
#' @param X Numeric \eqn{n \times p} matrix.
#' @param Y Numeric \eqn{n \times q} matrix with the same number of rows as
#'   \code{X}.
#' @param b Lower wrapping threshold. The reference value is 1.5.
#' @param c Upper wrapping threshold, with \eqn{0<b<c}. The reference value is
#'   4.
#' @param alpha Probability used for the chi-squared robust-distance cutoff.
#'   The reference value is 0.975.
#' @param q1,q2 Smooth wrapping constants. Their reference values are 1.540793
#'   and 0.8622731, respectively.
#' @param max_score_components Maximum number of canonical score pairs used for
#'   MCD reweighting. The reference value is 5.
#' @param lambda_cap Optional upper bound for each analytical shrinkage
#'   intensity. The default, 1, leaves the analytical estimate uncapped.
#' @param n_xi Deprecated compatibility argument. The consistency factor is now
#'   evaluated deterministically by numerical integration, so this argument is
#'   ignored.
#'
#' @details
#' Each column is centered by its median, scaled by its Gaussian-consistent MAD,
#' and transformed with a smooth redescending wrapping function. The Gaussian
#' consistency factor is evaluated deterministically rather than by Monte Carlo.
#'
#' An analytical Ledoit--Wolf-type shrinkage estimate regularizes the initial
#' joint correlation matrix. Robust distances are then computed from at most
#' \code{max_score_components} pairs of initial canonical scores. Binary weights
#' obtained from these distances are used to construct a second joint
#' correlation matrix, which is again shrinkage regularized before the final
#' canonical decomposition. If score-space MCD estimation is unavailable or too
#' few observations are retained, the initial regularized solution is returned.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{cor}: final canonical correlations.
#'   \item \code{a}, \code{b}: final canonical direction matrices for \code{X}
#'     and \code{Y}.
#'   \item \code{weights}: binary score-space observation weights.
#'   \item \code{retained}: number of observations with weight one.
#'   \item \code{reweighted_final}: whether the second weighted fit was used.
#'   \item \code{shrink_initial}, \code{shrink_final}: initial and final
#'     shrinkage intensities.
#'   \item \code{shrink_used}: backward-compatible alias of
#'     \code{shrink_final}.
#'   \item \code{xi}: deterministic Gaussian consistency factor.
#'   \item \code{wrapping}: wrapping constants used by the fit.
#'   \item \code{u1}, \code{v1}: first canonical score pair from the initial
#'     regularized solution, retained for backward compatibility and diagnostic
#'     plotting.
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 60
#' p <- 40
#' q <- 50
#' latent <- rnorm(n)
#' X <- outer(latent, seq(0.80, 0.30, length.out = p)) +
#'   matrix(rnorm(n * p), n, p)
#' Y <- outer(latent, seq(0.75, 0.25, length.out = q)) +
#'   matrix(rnorm(n * q), n, q)
#' fit <- AWRcca(X, Y)
#' fit$cor[1]
#' fit$retained
#' fit$shrink_final
#'
#' @importFrom MASS cov.mcd
#' @importFrom Matrix nearPD
#' @export
AWRcca <- function(X, Y, b = 1.5, c = 4, alpha = 0.975,
                   q1 = 1.540793, q2 = 0.8622731,
                   max_score_components = 5L,
                   lambda_cap = 1, n_xi = NULL) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("`X` and `Y` must be numeric matrices.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y)) {
    stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  }
  if (nrow(X) < 3L || ncol(X) < 1L || ncol(Y) < 1L) {
    stop("At least three observations and one variable per block are required.",
         call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("`X` and `Y` must contain only finite values.", call. = FALSE)
  }
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b) || b <= 0 ||
      !is.numeric(c) || length(c) != 1L || !is.finite(c) || c <= b) {
    stop("`b` and `c` must be finite scalars satisfying 0 < b < c.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a finite scalar in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(q1) || length(q1) != 1L || !is.finite(q1) || q1 <= 0 ||
      !is.numeric(q2) || length(q2) != 1L || !is.finite(q2) || q2 <= 0) {
    stop("`q1` and `q2` must be positive finite scalars.", call. = FALSE)
  }
  if (!is.numeric(max_score_components) || length(max_score_components) != 1L ||
      !is.finite(max_score_components) || max_score_components < 1) {
    stop("`max_score_components` must be a positive integer.", call. = FALSE)
  }
  max_score_components <- as.integer(max_score_components)
  if (!is.numeric(lambda_cap) || length(lambda_cap) != 1L ||
      !is.finite(lambda_cap) || lambda_cap <= 0 || lambda_cap > 1) {
    stop("`lambda_cap` must be a finite scalar in (0, 1].", call. = FALSE)
  }
  if (!is.null(n_xi)) {
    warning("`n_xi` is deprecated and ignored; `xi` is evaluated deterministically.",
            call. = FALSE)
  }

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  xi <- .awr_consistency_factor(b, c, q1, q2)
  Xw <- .awr_wrap_matrix(X, b, c, q1, q2, xi)
  Yw <- .awr_wrap_matrix(Y, b, c, q1, q2, xi)
  Z <- cbind(Xw, Yw)

  initial_cor <- .awr_shrink_correlation(Z, lambda_cap)
  initial <- .awr_cca_from_correlation(initial_cor$R, p, q)
  m <- min(max_score_components, ncol(initial$a), ncol(initial$b))
  scores <- cbind(
    Xw %*% initial$a[, seq_len(m), drop = FALSE],
    Yw %*% initial$b[, seq_len(m), drop = FALSE]
  )
  u1 <- as.numeric(scores[, 1L])
  v1 <- as.numeric(scores[, m + 1L])

  weights <- tryCatch({
    mcd <- MASS::cov.mcd(scores)
    rd2 <- stats::mahalanobis(scores, mcd$center, mcd$cov)
    as.integer(rd2 <= stats::qchisq(alpha, df = ncol(scores)))
  }, error = function(e) rep.int(1L, n))

  min_retained <- max(8L, ncol(scores) + 2L)
  use_final <- sum(weights) >= min_retained
  final <- initial
  final_lambda <- initial_cor$lambda

  if (use_final) {
    final_cor <- .awr_shrink_correlation(
      Z[weights == 1L, , drop = FALSE], lambda_cap
    )
    final <- .awr_cca_from_correlation(final_cor$R, p, q)
    final_lambda <- final_cor$lambda
  }

  list(
    cor = final$cor,
    a = final$a,
    b = final$b,
    weights = weights,
    retained = sum(weights),
    reweighted_final = use_final,
    shrink_initial = initial_cor$lambda,
    shrink_final = final_lambda,
    shrink_used = final_lambda,
    xi = xi,
    wrapping = c(b = b, c = c, q1 = q1, q2 = q2),
    u1 = u1,
    v1 = v1
  )
}

.awr_psi <- function(x, b, c, q1, q2) {
  ax <- abs(x)
  out <- x
  middle <- ax > b & ax < c
  out[middle] <- q1 * tanh(q2 * (c - ax[middle])) * sign(x[middle])
  out[ax >= c] <- 0
  out
}

.awr_consistency_factor <- function(b, c, q1, q2) {
  central <- stats::integrate(
    function(z) z^2 * stats::dnorm(z), 0, b,
    subdivisions = 1000L, rel.tol = 1e-11
  )$value
  shoulder <- stats::integrate(
    function(z) (q1 * tanh(q2 * (c - z)))^2 * stats::dnorm(z), b, c,
    subdivisions = 1000L, rel.tol = 1e-11
  )$value
  1 / sqrt(2 * (central + shoulder))
}

.awr_wrap_matrix <- function(M, b, c, q1, q2, xi) {
  out <- vapply(seq_len(ncol(M)), function(j) {
    column <- M[, j]
    center <- stats::median(column)
    scale <- stats::mad(column, constant = 1.482602218505602)
    if (!is.finite(scale) || scale < 1e-10) scale <- 1e-10
    xi * .awr_psi((column - center) / scale, b, c, q1, q2)
  }, numeric(nrow(M)))
  matrix(out, nrow = nrow(M), ncol = ncol(M))
}

.awr_standardize <- function(Z) {
  centered <- sweep(Z, 2L, colMeans(Z), "-")
  scales <- sqrt(colSums(centered^2) / nrow(centered))
  good <- is.finite(scales) & scales > 1e-12
  standardized <- matrix(0, nrow(Z), ncol(Z))
  if (any(good)) {
    standardized[, good] <- sweep(
      centered[, good, drop = FALSE], 2L, scales[good], "/"
    )
  }
  list(Z = standardized, good = good)
}

.awr_shrink_correlation <- function(Z, lambda_cap) {
  standardized <- .awr_standardize(Z)
  X <- standardized$Z
  n <- nrow(X)
  d <- ncol(X)
  S <- crossprod(X) / n
  diag(S)[standardized$good] <- 1
  S[!is.finite(S)] <- 0

  # Remove possible floating-point indefiniteness before shrinkage.
  eigenvalues <- tryCatch(
    eigen(S, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  if (any(!is.finite(eigenvalues)) || min(eigenvalues, na.rm = TRUE) < -1e-8) {
    S <- as.matrix(Matrix::nearPD(S, corr = TRUE)$mat)
  }

  target <- diag(d)
  phi_matrix <- crossprod(X^2) / n - S^2
  phi <- sum(pmax(phi_matrix, 0))
  denominator <- sum((S - target)^2)
  lambda <- if (!is.finite(denominator) || denominator <= 1e-14) {
    1
  } else {
    min(1, max(0, phi / (n * denominator)))
  }
  lambda <- min(lambda, lambda_cap)
  R <- (1 - lambda) * S + lambda * target
  diag(R) <- 1
  list(R = R, lambda = lambda)
}

.awr_inverse_sqrt <- function(S, floor_value = 1e-8) {
  decomposition <- eigen(S, symmetric = TRUE)
  values <- pmax(Re(decomposition$values), floor_value)
  decomposition$vectors %*%
    diag(1 / sqrt(values), length(values)) %*%
    t(decomposition$vectors)
}

.awr_cca_from_correlation <- function(R, p, q) {
  ix <- seq_len(p)
  iy <- p + seq_len(q)
  Rxx <- R[ix, ix, drop = FALSE]
  Ryy <- R[iy, iy, drop = FALSE]
  Rxy <- R[ix, iy, drop = FALSE]
  Wx <- .awr_inverse_sqrt(Rxx)
  Wy <- .awr_inverse_sqrt(Ryy)
  decomposition <- svd(Wx %*% Rxy %*% Wy)
  list(
    cor = pmin(1, pmax(0, Re(decomposition$d))),
    a = Wx %*% decomposition$u,
    b = Wy %*% decomposition$v
  )
}
