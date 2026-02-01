#' Adaptive Wrapped Robust Canonical Correlation Analysis (AWRcca)
#'
#' Implements an adaptive wrapped robust canonical correlation analysis procedure for
#' potentially contaminated high-dimensional data. The method applies columnwise robust
#' standardization and wrapping to mitigate cellwise outliers, uses a Fisher-consistency
#' correction, enforces positive semi-definiteness of the correlation matrix, applies
#' Ledoit--Wolf type shrinkage, and performs an MCD-based reweighting in the canonical
#' score space to downweight casewise outliers.
#'
#' @param X A numeric matrix of dimension \eqn{n \times p}.
#' @param Y A numeric matrix of dimension \eqn{n \times q}.
#' @param b Lower wrapping threshold (\eqn{b < c}). Default is 1.5.
#' @param c Upper wrapping threshold. Default is 4.
#' @param alpha Reweighting cutoff probability for chi-square threshold. Default is 0.975.
#' @param n_xi Monte Carlo sample size for the consistency correction. Default is 10000.
#' @param lambda_cap Upper bound for the shrinkage intensity. Default is 0.5.
#'
#' @details
#' The wrapping transformation is based on a smooth redescending function
#' \eqn{\psi_{b,c}} applied to robust z-scores (median/MAD). The shrinkage intensity is
#' estimated in a Ledoit--Wolf spirit and then capped by \code{lambda_cap} to avoid
#' overshrinkage.
#'
#' The function returns (i) canonical correlations, (ii) the shrinkage intensity used,
#' (iii) 0/1 reweighting indicators, and (iv) the first canonical score pair computed from
#' the initial solution (useful for diagnostic plots).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{cor}: vector of canonical correlations.
#'   \item \code{shrink_used}: shrinkage intensity used in the correlation regularization.
#'   \item \code{weights}: 0/1 weights from MCD-based reweighting in score space.
#'   \item \code{u1}: first canonical score for \code{X} (initial solution).
#'   \item \code{v1}: first canonical score for \code{Y} (initial solution).
#' }
#'
#' @examples
#' # Example: correlated blocks via a shared latent factor
#' set.seed(123)
#' n <- 50; p <- 30; q <- 20
#' u <- rnorm(n)
#' ax <- rnorm(p); ax <- ax / sqrt(sum(ax^2))
#' by <- rnorm(q); by <- by / sqrt(sum(by^2))
#' X <- 1.0 * u %*% t(ax) + matrix(rnorm(n*p), n, p)
#' Y <- 1.0 * u %*% t(by) + matrix(rnorm(n*q), n, q)
#' fit <- AWRcca(X, Y)
#' fit$cor[1]
#'
#' @importFrom Matrix nearPD
#' @importFrom MASS cov.mcd
#' @export
AWRcca <- function(X, Y, b = 1.5, c = 4.0,
                   alpha = 0.975, n_xi = 10000,
                   lambda_cap = 0.5) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  if (!is.numeric(X) || !is.numeric(Y)) stop("`X` and `Y` must be numeric matrices.")
  if (nrow(X) != nrow(Y)) stop("`X` and `Y` must have the same number of rows.")
  if (!(is.numeric(b) && is.numeric(c) && length(b) == 1L && length(c) == 1L && b < c)) {
    stop("`b` and `c` must be numeric scalars with b < c.")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0,1).")
  if (!is.numeric(n_xi) || n_xi < 1000) stop("`n_xi` should be >= 1000.")
  if (!is.numeric(lambda_cap) || lambda_cap < 0 || lambda_cap > 1) {
    stop("`lambda_cap` must be in [0,1].")
  }
  
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  
  # 1) Wrapping (columnwise)
  params <- .get_wrapping_params(b, c)
  
  wrap_fun <- function(Mat) {
    out <- apply(Mat, 2, function(col) {
      med <- stats::median(col)
      sig <- stats::mad(col)
      if (!is.finite(sig) || sig < 1e-9) sig <- 1e-6
      .psi_wrap((col - med) / sig, b, c, params$q1, params$q2)
    })
    # apply() returns a matrix with columns as variables; ensure dimension consistency
    out <- as.matrix(out)
    if (nrow(out) != n) out <- t(out)
    out
  }
  
  X_st <- wrap_fun(X)
  Y_st <- wrap_fun(Y)
  
  # 2) Consistency correction
  zz <- stats::rnorm(n_xi)
  vv <- stats::var(.psi_wrap(zz, b, c, params$q1, params$q2))
  xi <- 1 / sqrt(vv)
  if (!is.finite(xi) || xi <= 0) xi <- 1
  
  X_tilde <- xi * X_st
  Y_tilde <- xi * Y_st
  
  # 3) PSD & Shrinkage (correlation)
  Z <- cbind(X_tilde, Y_tilde)
  
  # handle zero-variance columns to avoid NA/Inf in cor()
  sdz  <- apply(Z, 2, stats::sd)
  good <- which(is.finite(sdz) & sdz > 1e-12)
  
  if (length(good) < 2) {
    R <- diag(ncol(Z))
  } else {
    Zg <- Z[, good, drop = FALSE]
    Rg <- stats::cor(Zg)
    Rg[!is.finite(Rg)] <- 0
    
    R <- diag(ncol(Z))
    R[good, good] <- Rg
  }
  
  # nearPD for numerical safety
  ev <- tryCatch(eigen(R, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
  if (any(!is.finite(ev)) || min(ev, na.rm = TRUE) < 1e-8) {
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    diag(R) <- 1
  }
  
  lam <- .estimate_lambda_lw(Z)
  lam <- min(max(lam, 0), lambda_cap)   # cap to avoid overshrinkage
  R_reg <- (1 - lam) * R + lam * diag(nrow(R))
  
  # 4) Initial CCA via SVD on block matrix
  Rxx <- R_reg[1:p, 1:p, drop = FALSE]
  Ryy <- R_reg[(p + 1):(p + q), (p + 1):(p + q), drop = FALSE]
  Rxy <- R_reg[1:p, (p + 1):(p + q), drop = FALSE]
  
  # Cholesky; fall back to identity if fails
  Lx <- tryCatch(t(chol(Rxx)), error = function(e) diag(p))
  Ly <- tryCatch(t(chol(Ryy)), error = function(e) diag(q))
  
  K  <- solve(Lx) %*% Rxy %*% t(solve(Ly))
  svd_res <- svd(K)
  
  a_init <- t(solve(Lx)) %*% svd_res$u
  b_init <- t(solve(Ly)) %*% svd_res$v
  
  # canonical scores for diagnostics (1st component) from initial solution
  u1_all <- as.numeric(X_tilde %*% a_init[, 1])
  v1_all <- as.numeric(Y_tilde %*% b_init[, 1])
  
  # 5) Reweighting (MCD on first few scores)
  n_comps <- min(5, ncol(a_init))
  scores <- cbind(
    X_tilde %*% a_init[, 1:n_comps, drop = FALSE],
    Y_tilde %*% b_init[, 1:n_comps, drop = FALSE]
  )
  
  w <- tryCatch({
    mcd <- MASS::cov.mcd(scores)
    rd2 <- stats::mahalanobis(scores, mcd$center, mcd$cov)
    cutoff <- stats::qchisq(alpha, df = ncol(scores))
    ifelse(rd2 <= cutoff, 1L, 0L)
  }, error = function(e) rep(1L, n))
  
  # 6) Final estimation (if feasible)
  if (sum(w) > (p + q + 2)) {
    res <- stats::cancor(X[w == 1, , drop = FALSE], Y[w == 1, , drop = FALSE])
    return(list(
      cor = res$cor,
      shrink_used = lam,
      weights = w,
      u1 = u1_all,
      v1 = v1_all
    ))
  }
  
  # fallback: use SVD canonical correlations
  list(
    cor = svd_res$d,
    shrink_used = lam,
    weights = w,
    u1 = u1_all,
    v1 = v1_all
  )
}

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' @keywords internal
.get_wrapping_params <- function(b, c) {
  f_opt <- function(q2) {
    A <- q2 * (c - b)
    val <- tanh(A)
    (b / val) * q2 * (1 - val^2) - 1
  }
  tryCatch({
    q2_sol <- stats::uniroot(f_opt, interval = c(1e-5, 10))$root
    q1_sol <- b / tanh(q2_sol * (c - b))
    list(q1 = q1_sol, q2 = q2_sol)
  }, error = function(e) list(q1 = b, q2 = 1))
}

#' @keywords internal
.psi_wrap <- function(x, b, c, q1, q2) {
  abs_x <- abs(x); sign_x <- sign(x); y <- x
  idx_mid <- (abs_x > b) & (abs_x < c)
  if (any(idx_mid)) {
    y[idx_mid] <- (q1 * tanh(q2 * (c - abs_x[idx_mid]))) * sign_x[idx_mid]
  }
  y[abs_x >= c] <- 0
  y
}

#' @keywords internal
.estimate_lambda_lw <- function(X) {
  n <- nrow(X); p <- ncol(X)
  if (n < 2) return(0.1)
  sample_cov <- stats::cov(X)
  prior <- mean(diag(sample_cov)) * diag(p)
  d <- mean((sample_cov - prior)^2)
  if (!is.finite(d) || d < 1e-10) d <- 1e-10
  r2 <- 1 / n * sum((crossprod(X^2) / n - sample_cov^2)^2)
  lam <- r2 / d
  lam <- max(0, min(1, lam))
  lam
}
