#' Adaptive Wrapped Robust Canonical Correlation Analysis (AWRcca)
#'
#' Implements the proposed Adaptive Wrapped Robust Canonical Correlation Analysis (AWRcca)
#' for high-dimensional data possibly contaminated by both cellwise and casewise outliers.
#' The method combines (i) columnwise robust standardization, (ii) wrapping transformation
#' for cellwise robustness, (iii) Fisher-consistency correction under normality,
#' (iv) PSD guarantee via nearest correlation matrix projection, (v) analytical Ledoit–Wolf
#' shrinkage regularization, and (vi) a reweighting step based on robust distances in the
#' canonical score space to downweight structural casewise outliers.
#'
#' @param X A numeric matrix of dimension \eqn{n \times p}.
#' @param Y A numeric matrix of dimension \eqn{n \times q}.
#' @param b Lower wrapping threshold (\eqn{b < c}). Default is 1.5.
#' @param c Upper wrapping threshold. Default is 4.
#' @param m_scores Number of canonical pairs used for reweighting in score space.
#'        Default is 5; internally truncated to \code{min(5, rank)}.
#' @param alpha Reweighting cutoff probability for chi-square threshold. Default is 0.975.
#' @param eps_psd Numerical tolerance for PSD checks. Default is 1e-10.
#' @param center_fun Robust centering function. Default is \code{stats::median}.
#' @param scale_fun Robust scale function. Default is MAD with consistency factor:
#'        \code{stats::mad(x, constant = 1.4826)}.
#' @param seed Optional integer seed for the Monte Carlo approximation used in the
#'        consistency correction factor \eqn{\xi}. Default is \code{NULL}.
#' @param xi_mc_n Monte Carlo size for computing \eqn{\xi}. Default is 200000.
#'        Increase for more precision.
#' @param return_all Logical; if TRUE, returns additional intermediate matrices (may be large).
#'        Default is FALSE.
#'        
#'       
#' @importFrom Matrix nearPD
#'
#' @details
#' **Overview**
#'
#' 1) **Robust standardization**: each column is centered by \code{median} and scaled by MAD.
#'
#' 2) **Wrapping**: robust z-scores are mapped by \eqn{\psi_{b,c}} (Raymaekers & Rousseeuw, 2021):
#' \deqn{
#' \psi_{b,c}(z) =
#' \begin{cases}
#' z, & |z| \le b, \\
#' q_1 \tanh(q_2 (c-|z|)) \,\mathrm{sign}(z), & b < |z| < c, \\
#' 0, & |z| \ge c,
#' \end{cases}
#' }
#' where \eqn{q_1 = b} and \eqn{q_2 = \tanh^{-1}(1)/ (c-b)} (a smooth choice ensuring continuity).
#'
#' 3) **Consistency correction**: a factor \eqn{\xi} is computed such that
#' \eqn{\mathrm{Var}(\xi \psi_{b,c}(Z)) = 1} for \eqn{Z \sim N(0,1)}.
#' Here \eqn{\xi = 1/\sqrt{\mathbb{E}[\psi_{b,c}(Z)^2]}} is estimated by Monte Carlo.
#'
#' 4) **Correlation & PSD**: a correlation matrix is formed from the transformed data.
#' If it is not PSD, it is projected via \code{Matrix::nearPD(..., corr = TRUE)}.
#'
#' 5) **Ledoit–Wolf shrinkage**: a data-driven shrinkage intensity \eqn{\lambda \in [0,1]}
#' is estimated for the (correlation) matrix with identity target and applied as
#' \eqn{R_\text{reg} = (1-\lambda)R + \lambda I}.
#'
#' 6) **CCA by SVD**: canonical correlations are obtained from the whitened cross-correlation block.
#'
#' 7) **Reweighting**: canonical scores are formed and robust distances are computed in the
#' concatenated score space using \code{robustbase::covMcd}. Observations with distances exceeding
#' the \eqn{\chi^2} cutoff are downweighted (0/1 weights).
#'
#' 8) **Final estimation**: a weighted correlation matrix is computed and steps 4–6 are repeated.
#'
#' @return
#' A list with class \code{"AWRcca"} containing:
#' \itemize{
#' \item \code{rho}: estimated canonical correlations.
#' \item \code{A}: canonical vectors for \code{X}.
#' \item \code{B}: canonical vectors for \code{Y}.
#' \item \code{weights}: 0/1 weights from the reweighting step.
#' \item \code{lambda}: shrinkage intensity used in the final estimation.
#' \item \code{lambda_init}: shrinkage intensity used in the initial estimation.
#' \item \code{xi}: consistency correction factor.
#' \item \code{b}, \code{c}, \code{alpha}, \code{m_scores}: used tuning parameters.
#' \item \code{n}, \code{p}, \code{q}: dimensions.
#' \item \code{runtime}: named vector with timing information (seconds).
#' \item \code{call}: matched call.
#' }
#' If \code{return_all = TRUE}, includes \code{R_init}, \code{R_final} and transformed data matrices.
#'
#' @examples
#' \dontrun{
#' # Example with random matrices
#' set.seed(1)
#' X <- matrix(rnorm(40*120), 40, 120)
#' Y <- matrix(rnorm(40*21),  40, 21)
#' fit <- AWRcca(X, Y)
#' fit
#' summary(fit)
#' }
#'
#' @export
AWRcca <- function(X, Y,
                   b = 1.5, c = 4,
                   m_scores = 5, alpha = 0.975,
                   eps_psd = 1e-10,
                   center_fun = stats::median,
                   scale_fun  = function(x) stats::mad(x, constant = 1.4826),
                   seed = NULL, xi_mc_n = 200000,
                   return_all = FALSE) {

  t0 <- proc.time()[3]

  # --- checks ---
  if (!is.matrix(X) || !is.numeric(X)) stop("`X` must be a numeric matrix.")
  if (!is.matrix(Y) || !is.numeric(Y)) stop("`Y` must be a numeric matrix.")
  if (nrow(X) != nrow(Y)) stop("`X` and `Y` must have the same number of rows (n).")
  if (!(is.numeric(b) && is.numeric(c) && length(b) == 1L && length(c) == 1L && b < c)) {
    stop("`b` and `c` must be numeric scalars with b < c.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required (for nearPD). Please install it.")
  }
  if (!requireNamespace("robustbase", quietly = TRUE)) {
    stop("Package 'robustbase' is required (for covMcd). Please install it.")
  }

  n <- nrow(X); p <- ncol(X); q <- ncol(Y)

  # --- Step 1: robust standardization + wrapping ---
  Zx <- .awrcca_wrap_block(X, b = b, c = c, center_fun = center_fun, scale_fun = scale_fun)
  Zy <- .awrcca_wrap_block(Y, b = b, c = c, center_fun = center_fun, scale_fun = scale_fun)

  # --- Step 2: consistency correction xi ---
  xi <- .awrcca_xi(b = b, c = c, mc_n = xi_mc_n, seed = seed)
  Zx <- xi * Zx
  Zy <- xi * Zy

  # --- Step 3-5: initial robust correlation -> PSD -> shrinkage ---
  R_init <- .awrcca_cor_safe(cbind(Zx, Zy))
  R_init <- .awrcca_make_psd(R_init, eps_psd = eps_psd)
  lambda_init <- .awrcca_ledoit_wolf_lambda(cbind(Zx, Zy), target = "identity")
  Rreg_init <- (1 - lambda_init) * R_init + lambda_init * diag(p + q)

  # --- Step 5: initial CCA via SVD on whitened cross-block ---
  cca_init <- .awrcca_svd_cca(Rreg_init, p = p, q = q)
  rho_init <- cca_init$rho
  A_init <- cca_init$A
  B_init <- cca_init$B

  t_init <- proc.time()[3]

  # --- Step 6: structural outlier detection (reweighting) in canonical score space ---
  # use first m_scores canonical pairs (or available)
  m_use <- min(m_scores, length(rho_init), 5L)
  if (m_use < 1L) {
    # fallback: no reweighting possible
    weights <- rep(1, n)
  } else {
    U <- Zx %*% A_init[, seq_len(m_use), drop = FALSE]
    V <- Zy %*% B_init[, seq_len(m_use), drop = FALSE]
    S_scores <- cbind(U, V) # n x (2m)
    weights <- .awrcca_reweight_scores(S_scores, alpha = alpha)
  }

  # --- Step 7: final estimation using weighted correlation + PSD + shrinkage + SVD ---
  R_final <- .awrcca_weighted_cor(cbind(Zx, Zy), w = weights)
  R_final <- .awrcca_make_psd(R_final, eps_psd = eps_psd)

  lambda_final <- .awrcca_ledoit_wolf_lambda(cbind(Zx, Zy), w = weights, target = "identity")
  Rreg_final <- (1 - lambda_final) * R_final + lambda_final * diag(p + q)

  cca_final <- .awrcca_svd_cca(Rreg_final, p = p, q = q)

  t_end <- proc.time()[3]

  out <- list(
    rho = cca_final$rho,
    A   = cca_final$A,
    B   = cca_final$B,
    weights = weights,
    lambda = lambda_final,
    lambda_init = lambda_init,
    xi = xi,
    b = b, c = c, alpha = alpha, m_scores = m_scores,
    n = n, p = p, q = q,
    runtime = c(
      total = t_end - t0,
      initial_fit = t_init - t0,
      reweight_and_final = t_end - t_init
    ),
    call = match.call()
  )

  if (isTRUE(return_all)) {
    out$R_init  <- Rreg_init
    out$R_final <- Rreg_final
    out$Zx <- Zx
    out$Zy <- Zy
  }

  class(out) <- "AWRcca"
  out
}

#' @export
print.AWRcca <- function(x, ...) {
  cat("AWRcca fit\n")
  cat("  n =", x$n, " p =", x$p, " q =", x$q, "\n")
  cat("  shrinkage lambda (final) =", formatC(x$lambda, digits = 4, format = "f"), "\n")
  cat("  inliers =", sum(x$weights), "/", length(x$weights), "\n")
  k <- min(5, length(x$rho))
  cat("  first canonical correlations:\n")
  cat("   ", paste(formatC(x$rho[seq_len(k)], digits = 4, format = "f"), collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.AWRcca <- function(object, ...) {
  s <- list(
    call = object$call,
    dims = c(n = object$n, p = object$p, q = object$q),
    rho = object$rho,
    lambda = object$lambda,
    lambda_init = object$lambda_init,
    xi = object$xi,
    inliers = sum(object$weights),
    outliers = length(object$weights) - sum(object$weights),
    runtime = object$runtime
  )
  class(s) <- "summary.AWRcca"
  s
}

#' @export
print.summary.AWRcca <- function(x, ...) {
  cat("Summary of AWRcca fit\n")
  cat("  n =", x$dims["n"], " p =", x$dims["p"], " q =", x$dims["q"], "\n")
  cat("  xi =", formatC(x$xi, digits = 5, format = "f"),
      " lambda_init =", formatC(x$lambda_init, digits = 4, format = "f"),
      " lambda_final =", formatC(x$lambda, digits = 4, format = "f"), "\n")
  cat("  inliers =", x$inliers, " outliers =", x$outliers, "\n")
  k <- min(10, length(x$rho))
  cat("  canonical correlations (first", k, "):\n")
  cat("   ", paste(formatC(x$rho[seq_len(k)], digits = 4, format = "f"), collapse = ", "), "\n")
  cat("  runtime (s):\n")
  cat("    total =", formatC(x$runtime["total"], digits = 3, format = "f"),
      " initial =", formatC(x$runtime["initial_fit"], digits = 3, format = "f"),
      " final =", formatC(x$runtime["reweight_and_final"], digits = 3, format = "f"), "\n")
  invisible(x)
}

# =============================================================================
# Internal helpers
# =============================================================================

#' Wrapping function \eqn{\psi_{b,c}}
#' @keywords internal
.awrcca_psi <- function(z, b, c) {
  az <- abs(z)
  s  <- sign(z)

  # Smooth transition parameters (simple, continuous choice)
  # q1 = b, q2 chosen so that tanh(q2*(c-b)) ~ 1
  q1 <- b
  q2 <- atanh(0.999999) / (c - b)

  out <- z
  mid <- (az > b) & (az < c)
  out[mid] <- q1 * tanh(q2 * (c - az[mid])) * s[mid]
  out[az >= c] <- 0
  out
}

#' Wrap a data block columnwise after robust standardization
#' @keywords internal
.awrcca_wrap_block <- function(M, b, c, center_fun, scale_fun) {
  n <- nrow(M); p <- ncol(M)
  out <- matrix(0, n, p)
  for (j in seq_len(p)) {
    x <- M[, j]
    mu <- center_fun(x)
    sc <- scale_fun(x)
    if (!is.finite(sc) || sc <= 0) sc <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sc) || sc <= 0) sc <- 1
    z <- (x - mu) / sc
    out[, j] <- .awrcca_psi(z, b = b, c = c)
  }
  colnames(out) <- colnames(M)
  out
}

#' Fisher-consistency correction factor xi under N(0,1)
#' xi = 1 / sqrt(E[psi(Z)^2])
#' @keywords internal
.awrcca_xi <- function(b, c, mc_n = 200000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  z <- stats::rnorm(mc_n)
  w <- .awrcca_psi(z, b = b, c = c)
  e2 <- mean(w^2)
  if (!is.finite(e2) || e2 <= 0) return(1)
  1 / sqrt(e2)
}

#' Safe correlation computation
#' @keywords internal
.awrcca_cor_safe <- function(Z) {
  # Use covariance then standardize to correlation for stability
  S <- stats::cov(Z)
  d <- sqrt(diag(S))
  d[d == 0] <- 1
  R <- S / (d %o% d)
  diag(R) <- 1
  R
}

#' Ensure PSD correlation matrix (nearPD if needed)
#' @keywords internal
.awrcca_make_psd <- function(R, eps_psd = 1e-10) {
  ev <- tryCatch(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values),
                 error = function(e) NA_real_)
  if (!is.finite(ev) || ev < eps_psd) {
    Rnpd <- Matrix::nearPD(R, corr = TRUE)$mat
    Rnpd <- as.matrix(Rnpd)
    diag(Rnpd) <- 1
    return(Rnpd)
  }
  R
}

#' Ledoit–Wolf shrinkage intensity for correlation/cov with identity target
#'
#' We estimate lambda = min(max(beta_hat / delta_hat, 0), 1) where
#' beta_hat estimates variance of sample covariance entries and
#' delta_hat = ||S - I||_F^2.
#'
#' @keywords internal
.awrcca_ledoit_wolf_lambda <- function(Z, w = NULL, target = c("identity")) {
  target <- match.arg(target)

  n <- nrow(Z)
  if (is.null(w)) {
    w <- rep(1, n)
  } else {
    w <- as.numeric(w)
    if (length(w) != n) stop("weights length mismatch.")
  }
  # normalize weights to sum to n_eff for stability
  w <- w / mean(w)
  n_eff <- sum(w > 0)

  # Weighted centering (Z already approximately centered via robust standardization,
  # but we center again for LW estimation)
  mu <- colSums(Z * w) / sum(w)
  Zc <- sweep(Z, 2, mu, "-")

  # Weighted covariance (unbiased-ish scaling)
  # Use denom = sum(w) - 1 to mimic sample covariance
  denom <- max(sum(w) - 1, 1)
  S <- crossprod(Zc * sqrt(w)) / denom

  # Convert to correlation-like scale if diagonal deviates
  d <- sqrt(diag(S)); d[d == 0] <- 1
  R <- S / (d %o% d)
  diag(R) <- 1

  # delta_hat = ||R - I||_F^2
  delta_hat <- sum((R - diag(ncol(R)))^2)

  if (!is.finite(delta_hat) || delta_hat <= 0) return(0)

  # beta_hat: average Frobenius norm of (x_i x_i' - S) squared
  # For correlation target, use standardized rows
  Zstd <- sweep(Zc, 2, sqrt(diag(S)), "/")
  Zstd[!is.finite(Zstd)] <- 0

  # compute beta_hat efficiently
  beta_sum <- 0
  for (i in seq_len(n)) {
    if (w[i] <= 0) next
    xi <- Zstd[i, , drop = FALSE]
    Si <- crossprod(xi) # outer product
    beta_sum <- beta_sum + w[i] * sum((Si - R)^2)
  }
  beta_hat <- beta_sum / (sum(w)^2) * n_eff  # scaled to be roughly comparable

  lam <- beta_hat / delta_hat
  lam <- max(min(lam, 1), 0)
  lam
}

#' Solve CCA via SVD from a regularized correlation matrix
#' @keywords internal
.awrcca_svd_cca <- function(Rreg, p, q) {
  Rxx <- Rreg[seq_len(p), seq_len(p), drop = FALSE]
  Ryy <- Rreg[p + seq_len(q), p + seq_len(q), drop = FALSE]
  Rxy <- Rreg[seq_len(p), p + seq_len(q), drop = FALSE]

  # Cholesky; fallback to nearPD if needed
  Lx <- tryCatch(chol(Rxx), error = function(e) NULL)
  Ly <- tryCatch(chol(Ryy), error = function(e) NULL)

  if (is.null(Lx)) {
    Rxx2 <- .awrcca_make_psd(Rxx)
    Lx <- chol(Rxx2)
  }
  if (is.null(Ly)) {
    Ryy2 <- .awrcca_make_psd(Ryy)
    Ly <- chol(Ryy2)
  }

  # Whitening: K = Lx^{-1} Rxy Ly^{-T}
  Lx_inv <- backsolve(Lx, diag(p), upper.tri = TRUE)
  Ly_inv <- backsolve(Ly, diag(q), upper.tri = TRUE)

  K <- Lx_inv %*% Rxy %*% t(Ly_inv)

  s <- svd(K)
  rho <- pmin(pmax(s$d, 0), 1)

  A <- t(Lx_inv) %*% s$u
  B <- t(Ly_inv) %*% s$v

  colnames(A) <- paste0("can", seq_len(ncol(A)))
  colnames(B) <- paste0("can", seq_len(ncol(B)))

  list(rho = rho, A = A, B = B)
}

#' Robust reweighting in score space using MCD distances
#' @keywords internal
.awrcca_reweight_scores <- function(S, alpha = 0.975) {
  n <- nrow(S)
  d <- ncol(S)
  if (d < 1) return(rep(1, n))

  # robustbase::covMcd works well for small d (here d = 2m <= 10)
  fit <- robustbase::covMcd(S)
  mu <- fit$center
  Sigma <- fit$cov
  # numerical safety
  Sigma <- .awrcca_make_psd(Sigma, eps_psd = 1e-12)

  rd2 <- stats::mahalanobis(S, center = mu, cov = Sigma)
  cutoff <- stats::qchisq(alpha, df = d)
  as.integer(rd2 <= cutoff)
}

#' Weighted correlation matrix with 0/1 weights
#' @keywords internal
.awrcca_weighted_cor <- function(Z, w) {
  n <- nrow(Z)
  w <- as.numeric(w)
  if (length(w) != n) stop("weights length mismatch.")
  w[w < 0] <- 0
  if (sum(w) <= 1) {
    # degenerate: fallback to unweighted
    return(.awrcca_cor_safe(Z))
  }
  # normalize weights to mean 1
  w <- w / mean(w)

  mu <- colSums(Z * w) / sum(w)
  Zc <- sweep(Z, 2, mu, "-")
  denom <- max(sum(w) - 1, 1)
  S <- crossprod(Zc * sqrt(w)) / denom

  d <- sqrt(diag(S))
  d[d == 0] <- 1
  R <- S / (d %o% d)
  diag(R) <- 1
  R
}
