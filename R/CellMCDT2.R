#' @title Cellwise Robust Permutation Hotelling T^2 Test for Two Independent Samples
#'
#' @description
#' Performs the cellMCD-based robust permutation Hotelling T^2 test for comparing
#' the mean vectors of two independent multivariate samples in the presence of
#' cellwise outliers.
#'
#' @details
#' The \\code{CellMCDT2} function implements the robust two-sample Hotelling
#' T^2 procedure based on the cellwise minimum covariance determinant
#' (cellMCD) estimator. For each group, the cellMCD location vector and
#' covariance matrix are computed. These estimates are then used to construct
#' a robust pooled covariance matrix and the observed robust Hotelling T^2
#' statistic. Since the finite-sample null distribution of the statistic is
#' unknown, the p-value is obtained by randomly permuting the group labels.
#'
#' The test is designed for settings in which individual cells of the data
#' matrix may be contaminated. It can also be used under clean data, where it
#' typically provides competitive performance while protecting the inference
#' against cellwise outliers.
#'
#' @param X1 A numeric matrix or data frame containing the observations from
#' the first group. Rows are observations and columns are variables.
#' @param X2 A numeric matrix or data frame containing the observations from
#' the second group. Rows are observations and columns are variables.
#' @param B The number of random permutations used to approximate the null
#' distribution. The default is \\code{B = 100}.
#' @param alpha The cellMCD trimming parameter. In each column, at least
#' \\code{n * alpha} cells must remain unflagged. The default is
#' \\code{alpha = 0.75}.
#' @param quant Cutoff probability used by \\code{cellWise::cellMCD} to flag
#' outlying cells. The default is \\code{quant = 0.99}.
#' @param crit Convergence tolerance used by \\code{cellWise::cellMCD}. The
#' default is \\code{crit = 1e-4}.
#' @param noCits Maximum number of C-steps used by \\code{cellWise::cellMCD}.
#' The default is \\code{noCits = 100}.
#' @param lmin Lower bound on the eigenvalues of the estimated covariance
#' matrix on the standardized data. The default is \\code{lmin = 1e-4}.
#' @param seed Optional integer used to set the random seed for the permutation
#' procedure. The default is \\code{NULL}.
#' @param correction Logical. If \\code{TRUE}, the permutation p-value is
#' calculated as \\code{(sum(T_perm >= T_obs) + 1) / (B + 1)}. If \\code{FALSE},
#' it is calculated as \\code{sum(T_perm >= T_obs) / B}. The default is
#' \\code{TRUE}.
#' @param verbose Logical. If \\code{TRUE}, progress information is printed
#' during the permutation procedure. The default is \\code{FALSE}.
#'
#' @return A list of class \\code{MVTests} with the following elements:
#' \\item{T2}{The observed cellMCD-based robust Hotelling T^2 statistic.}
#' \\item{p.value}{The permutation-based p-value.}
#' \\item{Permutations_T2}{The test statistic values obtained from permutations.}
#' \\item{B}{The number of permutations.}
#' \\item{alpha}{The cellMCD trimming parameter.}
#' \\item{quant}{The cutoff probability used by cellMCD.}
#' \\item{n1}{The sample size of the first group.}
#' \\item{n2}{The sample size of the second group.}
#' \\item{p}{The number of variables.}
#' \\item{group.centers}{The cellMCD location estimates for the two groups.}
#' \\item{group.covariances}{The cellMCD covariance estimates for the two groups.}
#' \\item{pooled.covariance}{The robust pooled covariance matrix.}
#' \\item{Test}{The name of the test.}
#'
#' @references
#' Raymaekers, J. and Rousseeuw, P. J. (2024). The cellwise minimum covariance
#' determinant estimator. \\emph{Journal of the American Statistical Association},
#' 119(548), 2610--2621.
#'
#' Bulut, H. and Esmeray, M. A cellwise robust Hotelling test for two-sample comparisons (Unpublished).
#'
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#'
#' @export
#'
#' @examples
#' if (requireNamespace("mvtnorm", quietly = TRUE) &&
#'     requireNamespace("cellWise", quietly = TRUE)) {
#'   set.seed(123)
#'   x1 <- mvtnorm::rmvnorm(n = 30, mean = rep(0, 5), sigma = diag(5))
#'   x2 <- mvtnorm::rmvnorm(n = 30, mean = rep(0, 5), sigma = diag(5))
#'   fit <- CellMCDT2(X1 = x1, X2 = x2, B = 19, seed = 123)
#'   fit$p.value
#' }
CellMCDT2 <- function(X1, X2, B = 100, alpha = 0.75, quant = 0.99,
                      crit = 1e-4, noCits = 100, lmin = 1e-4,
                      seed = NULL, correction = TRUE, verbose = FALSE) {

  if (!requireNamespace("cellWise", quietly = TRUE)) {
    stop("Package 'cellWise' is required for CellMCDT2. Please install it.",
         call. = FALSE)
  }

  if (!is.matrix(X1) && !is.data.frame(X1)) {
    stop("'X1' must be a matrix or data frame.", call. = FALSE)
  }
  if (!is.matrix(X2) && !is.data.frame(X2)) {
    stop("'X2' must be a matrix or data frame.", call. = FALSE)
  }

  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  if (!is.numeric(X1) || !is.numeric(X2)) {
    stop("'X1' and 'X2' must contain numeric values only.", call. = FALSE)
  }
  if (anyNA(X1) || anyNA(X2)) {
    stop("Missing values are not allowed in 'X1' or 'X2'. Please remove or impute them.",
         call. = FALSE)
  }
  if (ncol(X1) != ncol(X2)) {
    stop("'X1' and 'X2' must have the same number of columns.", call. = FALSE)
  }

  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p <- ncol(X1)

  if (n1 < 3 || n2 < 3) {
    stop("Each group must contain at least three observations.", call. = FALSE)
  }
  if (!is.numeric(B) || B < 1) {
    stop("'B' must be a positive integer.", call. = FALSE)
  }
  B <- as.integer(B)
  if (!is.numeric(alpha) || alpha <= 0.5 || alpha > 1) {
    stop("'alpha' must be a numeric value in the interval (0.5, 1].", call. = FALSE)
  }
  if (!is.numeric(quant) || quant <= 0 || quant >= 1) {
    stop("'quant' must be a numeric value in the interval (0, 1).", call. = FALSE)
  }
  if (!is.logical(correction) || length(correction) != 1) {
    stop("'correction' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  observed <- .CellMCDT2_stat(X1 = X1, X2 = X2, alpha = alpha,
                              quant = quant, crit = crit,
                              noCits = noCits, lmin = lmin)
  T2_obs <- observed$T2
  perm_stats <- numeric(B)
  pooled <- rbind(X1, X2)
  n <- n1 + n2

  for (b in seq_len(B)) {
    ind1 <- sample.int(n = n, size = n1, replace = FALSE)
    X1b <- pooled[ind1, , drop = FALSE]
    X2b <- pooled[-ind1, , drop = FALSE]

    perm_stats[b] <- .CellMCDT2_stat(X1 = X1b, X2 = X2b,
                                     alpha = alpha, quant = quant,
                                     crit = crit, noCits = noCits,
                                     lmin = lmin)$T2

    if (verbose && b %% 10 == 0) {
      message("Permutation ", b, " of ", B, " completed.")
    }
  }

  if (correction) {
    pval <- (sum(perm_stats >= T2_obs) + 1) / (B + 1)
  } else {
    pval <- sum(perm_stats >= T2_obs) / B
  }

  results <- list(
    T2 = T2_obs,
    p.value = pval,
    Permutations_T2 = perm_stats,
    B = B,
    alpha = alpha,
    quant = quant,
    crit = crit,
    noCits = noCits,
    lmin = lmin,
    correction = correction,
    n1 = n1,
    n2 = n2,
    p = p,
    group.centers = list(Group1 = observed$mu1, Group2 = observed$mu2),
    group.covariances = list(Group1 = observed$S1, Group2 = observed$S2),
    pooled.covariance = observed$C,
    Test = "CellMCDT2"
  )

  class(results) <- c("MVTests", "list")
  return(results)
}

.CellMCDT2_stat <- function(X1, X2, alpha, quant, crit, noCits, lmin) {

  fit1 <- cellWise::cellMCD(X1, alpha = alpha, quant = quant,
                            crit = crit, noCits = noCits, lmin = lmin)
  fit2 <- cellWise::cellMCD(X2, alpha = alpha, quant = quant,
                            crit = crit, noCits = noCits, lmin = lmin)

  mu1 <- as.numeric(fit1$mu)
  mu2 <- as.numeric(fit2$mu)
  S1 <- as.matrix(fit1$S)
  S2 <- as.matrix(fit2$S)

  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n <- n1 + n2

  C <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  diff <- matrix(mu1 - mu2, ncol = 1)

  invC <- tryCatch(
    solve(C),
    error = function(e) {
      as.matrix(solve(as.matrix(Matrix::nearPD(C)$mat)))
    }
  )

  T2 <- as.numeric((n1 * n2 / n) * t(diff) %*% invC %*% diff)

  return(list(T2 = T2, mu1 = mu1, mu2 = mu2, S1 = S1, S2 = S2, C = C))
}
