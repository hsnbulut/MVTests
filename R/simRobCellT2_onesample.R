#' Monte Carlo Calibration for the One-Sample cellMCD Hotelling Test
#'
#' @description
#' Estimates d and q for an empirical scaled-F approximation to the one-sample
#' cellMCD Hotelling T^2 statistic under a Gaussian null model with fixed mean mu0.
#'
#' @param X Optional finite numeric matrix or data frame of observed data.
#'   If supplied, its row and column counts determine n and p. Missing values
#'   are rejected; perform complete-case handling before calling this function
#'   and use the same retained observations in the test.
#' @param mu0 Required finite numeric null mean vector of length p.
#' @param Sigma Optional finite symmetric positive-definite p-by-p covariance
#'   matrix. If supplied it takes precedence over covariance estimation from X.
#'   Otherwise X is required and the bias-corrected scatter S from
#'   cellWise::cellMCD(X) is used.
#' @param n Optional integer sample size greater than 1. Required without X;
#'   when supplied with X it must equal nrow(X).
#' @param p Optional integer dimension greater than 1. Required without X;
#'   when supplied with X it must equal ncol(X).
#' @param nper Positive integer number of Monte Carlo replications; default 3000.
#'
#' @details
#' When Sigma is omitted, cellMCD is applied to the raw X with its default
#' settings and a free location. This differs from the fixed-center plug-in
#' design of the earlier data-based interface. The observed fitted location is
#' never used as the simulation mean: every Monte Carlo sample is generated from
#' N_p(mu0, Sigma), and its cellMCD location estimate is compared with the fixed
#' mu0 in the statistic \eqn{n (location - mu0)^T S^{-1} (location - mu0)}.
#'
#' All cellMCD fits use the installed API's default settings, including
#' alpha = 0.75, quant = 0.99 and crit = 1e-04. Match those settings and the
#' effective n and p when applying RobCellT2_onesample(). Actual settings are
#' recorded in the result. Use set.seed() before calling for reproducibility.
#'
#' The d and q values match the Monte Carlo mean and variance to d * F(p, q).
#' They define an empirical approximation, not a theoretically exact or
#' universally covariance-invariant distribution. A plug-in covariance is an
#' estimate, not guaranteed to recover the population covariance exactly.
#' At least 10 positive finite statistics are required. Failed fits or invalid
#' statistics are recorded, and returned moments use successful replications.
#' Unstable moments, including matches that cannot yield finite q > 4, produce
#' an error. Inspect failure counts and consider increasing nper.
#'
#' The pre-2.3.6 n/p/nrep/seed/data interface is available as
#' simRobCellT2_onesample_legacy(), with a deprecation warning and its original
#' semantics. This function's six-argument signature does not accept nrep, seed,
#' data or additional cellMCD arguments. Migration requires explicit mu0; use
#' nper instead of nrep, X instead of data, and set.seed() instead of seed.
#'
#' @return A list containing d, q, mean.T2, var.T2, Sigma, mu0, n, p, nper,
#'   n.success and n.failed. It also contains failure.messages, covariance.source,
#'   mode, cellMCD.settings, covariance.settings and cellWise.version. The nrep
#'   element is retained as an alias for nper for consumers of earlier results.
#' @seealso \code{\link{RobCellT2_onesample}},
#'   \code{\link{simRobCellT2_onesample_legacy}}
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#' @examples
#' if (requireNamespace("MASS", quietly = TRUE) &&
#'     requireNamespace("cellWise", quietly = TRUE)) {
#'   mu0 <- c(1, -1)
#'   Sigma <- matrix(c(2, 0.4, 0.4, 1), 2)
#'   set.seed(123)
#'   X <- MASS::mvrnorm(30, mu0, Sigma)
#'   set.seed(1)
#'   constants <- simRobCellT2_onesample(X = X, mu0 = mu0, nper = 50)
#'   set.seed(1)
#'   simRobCellT2_onesample(mu0 = mu0, Sigma = Sigma, n = 30, p = 2, nper = 50)
#' }
#' @export
simRobCellT2_onesample <- function(X = NULL, mu0, Sigma = NULL,
                                  n = NULL, p = NULL, nper = 3000) {
  if (missing(mu0) || is.null(mu0)) {
    stop("mu0 is required: supply the fixed null mean vector.", call. = FALSE)
  }
  if (!is.null(X)) {
    if (!is.matrix(X) && !is.data.frame(X)) {
      stop("X must be a finite numeric matrix or data frame. Old positional n,p calls must use simRobCellT2_onesample_legacy().",
           call. = FALSE)
    }
    X <- as.matrix(X)
    if (!is.numeric(X) || any(!is.finite(X))) {
      stop("X must be a finite numeric matrix or data frame without missing values.",
           call. = FALSE)
    }
    if (!is.null(n) && (!.cellT2_integer(n, 2) || n != nrow(X))) {
      stop("n must be an integer greater than 1 and equal nrow(X).", call. = FALSE)
    }
    if (!is.null(p) && (!.cellT2_integer(p, 2) || p != ncol(X))) {
      stop("p must be an integer greater than 1 and equal ncol(X).", call. = FALSE)
    }
    n <- nrow(X)
    p <- ncol(X)
  } else if (is.null(n) || is.null(p) || is.null(Sigma)) {
    stop("Without X, n, p and Sigma must all be supplied.", call. = FALSE)
  }
  if (!.cellT2_integer(n, 2)) stop("n must be a finite integer greater than 1.", call. = FALSE)
  if (!.cellT2_integer(p, 2)) stop("p must be a finite integer greater than 1.", call. = FALSE)
  if (!.cellT2_integer(nper, 1)) stop("nper must be a positive finite integer.", call. = FALSE)
  n <- as.integer(n); p <- as.integer(p); nper <- as.integer(nper)
  if (!is.numeric(mu0) || !is.null(dim(mu0)) || length(mu0) != p || any(!is.finite(mu0))) {
    stop("mu0 must be a finite numeric vector of length p.", call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) stop("The 'MASS' package is required.", call. = FALSE)
  if (!requireNamespace("cellWise", quietly = TRUE)) stop("The 'cellWise' package is required.", call. = FALSE)
  estimated <- is.null(Sigma)
  if (estimated) {
    fit <- NULL
    invisible(utils::capture.output({
      fit <- tryCatch(cellWise::cellMCD(X), error = function(e) {
        stop(paste("cellMCD covariance estimation on raw X failed:", conditionMessage(e)),
             call. = FALSE)
      })
    }))
    if (is.null(fit) || is.null(fit$S)) {
      stop("cellMCD covariance estimation did not return the bias-corrected scatter S.",
           call. = FALSE)
    }
    Sigma <- .cellT2_covariance(fit$S, p, "cellMCD covariance S")
  } else Sigma <- .cellT2_covariance(Sigma, p)
  # Share the existing simulation kernel: free-location cellMCD fits and the
  # identical n * delta' solve(S) delta statistic, centered on the fixed mu0.
  result <- tryCatch(
    .simRobCellT2_onesample_legacy_impl(n = n, p = p, nrep = nper,
                                     Sigma = Sigma, mu0 = as.numeric(mu0)),
    error = function(e) stop(gsub("nrep", "nper", conditionMessage(e), fixed = TRUE),
                             call. = FALSE)
  )
  result$nper <- nper
  result$mode <- if (estimated) "data" else "explicit"
  result$covariance.source <- if (estimated) "cellMCD S from raw X" else "supplied Sigma"
  result$covariance.settings <- if (estimated) result$cellMCD.settings else NULL
  result
}

.cellT2_integer <- function(x, minimum) {
  is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= minimum && x <= .Machine$integer.max && x == floor(x)
}

.cellT2_covariance <- function(Sigma, p, label = "Sigma") {
  if (!is.matrix(Sigma) || !is.numeric(Sigma) ||
      !identical(dim(Sigma), c(p, p)) || any(!is.finite(Sigma))) {
    stop(paste(label, "must be a finite numeric p-by-p matrix."), call. = FALSE)
  }
  if (!isSymmetric(Sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
    stop(paste(label, "must be symmetric within numerical tolerance."), call. = FALSE)
  }
  Sigma <- Sigma / 2 + t(Sigma) / 2
  if (is.null(tryCatch(chol(Sigma), error = function(e) NULL))) {
    stop(paste(label, "must be positive definite."), call. = FALSE)
  }
  Sigma
}
