#' @title Weighted MRCD-Based Robust MANOVA Test for High-Dimensional Data
#'
#' @description
#' Performs a weighted minimum regularized covariance determinant (MRCD)-based
#' robust one-way MANOVA test for high-dimensional data.
#'
#' @details
#' The \code{RobHDMANOVA} function tests the equality of multivariate group
#' location vectors in one-way MANOVA settings, particularly when the number of
#' variables is large relative to the sample size and the data may contain
#' outlying observations.
#'
#' The procedure first computes groupwise MRCD location estimates. Then, a pooled
#' MRCD covariance matrix is obtained from group-centered observations. Robust
#' distances are calculated using this pooled covariance matrix, and binary
#' weights are assigned according to a robust distance cutoff. Reweighted group
#' means are then used to construct a robust between-group scatter matrix. A
#' robust Wilks-type statistic is computed as
#' \deqn{
#'   \Lambda_R = \frac{|W_R|}{|W_R + B_R|},
#' }
#' where \eqn{W_R} is the pooled MRCD covariance matrix and \eqn{B_R} is the
#' robust between-group scatter matrix. The test statistic is
#' \deqn{
#'   T_R = -\log(\Lambda_R).
#' }
#' Since the finite-sample null distribution is unknown, the p-value is obtained
#' using a permutation procedure.
#'
#' @param x A numeric data matrix or data frame. Rows represent observations and
#' columns represent variables.
#' @param group A grouping vector indicating the group membership of each
#' observation. It will be internally converted to a factor.
#' @param N The number of permutations used to approximate the null distribution.
#' The default is \code{N = 100}.
#' @param alpha Numeric parameter controlling the size of the subsets over which
#' the MRCD determinant is minimized. Allowed values are between 0.5 and 1.
#' The default is \code{alpha = 0.75}.
#' @param tau Cutoff probability used in the robust distance-based reweighting
#' step. The default is \code{tau = 0.975}.
#' @param cutoff The cutoff rule for robust distances. Options are
#' \code{"normal"} and \code{"chisq"}. The default is \code{"normal"}.
#' @param seed An optional integer used to set the random seed for the permutation
#' procedure. The default is \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, progress information is printed during
#' the permutation procedure. The default is \code{FALSE}.
#'
#' @return A list of class \code{MVTests} with the following elements:
#' \item{Lambda}{The robust Wilks' Lambda value.}
#' \item{TR}{The observed robust MANOVA test statistic.}
#' \item{p.value}{The permutation-based p-value.}
#' \item{Permutations_TR}{The test statistic values obtained from permutations.}
#' \item{alpha}{The trimming parameter used in MRCD estimation.}
#' \item{tau}{The cutoff probability used for robust distance-based reweighting.}
#' \item{cutoff}{The cutoff rule used for robust distances.}
#' \item{group.centers}{The reweighted robust group centers.}
#' \item{weights}{The binary robust weights for observations in each group.}
#' \item{Test}{The name of the test.}
#'
#' @references
#' Boudt, K., Rousseeuw, P. J., Vanduffel, S., and Verdonck, T. (2020).
#' The minimum regularized covariance determinant estimator.
#' \emph{Statistics and Computing}, 30, 113--128.
#'
#' Todorov, V. and Filzmoser, P. (2010).
#' Robust statistic for the one-way MANOVA.
#' \emph{Computational Statistics and Data Analysis}, 54, 37--48.
#'
#' Bulut, H. (2020).
#' Mahalanobis distance based on minimum regularized covariance determinant
#' estimators for high dimensional data.
#' \emph{Communications in Statistics - Theory and Methods}, 49, 5897--5907.
#'
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#'
#' @export
#'
#' @examples
#' if (requireNamespace("rrcov", quietly = TRUE) &&
#'     requireNamespace("mvtnorm", quietly = TRUE)) {
#'
#'   set.seed(123)
#'   x1 <- mvtnorm::rmvnorm(n = 10, mean = rep(0, 20), sigma = diag(20))
#'   x2 <- mvtnorm::rmvnorm(n = 10, mean = rep(0, 20), sigma = diag(20))
#'   x3 <- mvtnorm::rmvnorm(n = 10, mean = rep(0, 20), sigma = diag(20))
#'
#'   x <- rbind(x1, x2, x3)
#'   group <- c(rep(1, 10), rep(2, 10), rep(3, 10))
#'
#'   RobHDMANOVA(x = x, group = group, N = 19, alpha = 0.75,
#'               tau = 0.975, seed = 123)
#' }
RobHDMANOVA <- function(x, group, N = 100, alpha = 0.75, tau = 0.975,
                        cutoff = c("normal", "chisq"),
                        seed = NULL, verbose = FALSE) {
  
  if (!requireNamespace("rrcov", quietly = TRUE)) {
    stop("Package 'rrcov' is required for RobHDMANOVA. Please install it.",
         call. = FALSE)
  }
  
  cutoff <- match.arg(cutoff)
  
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.", call. = FALSE)
  }
  
  x <- as.matrix(x)
  
  if (!is.numeric(x)) {
    stop("'x' must contain numeric values only.", call. = FALSE)
  }
  
  if (anyNA(x)) {
    stop("'x' contains missing values. Please remove or impute missing values.",
         call. = FALSE)
  }
  
  group <- as.factor(group)
  
  if (nrow(x) != length(group)) {
    stop("The number of rows of 'x' must be equal to the length of 'group'.",
         call. = FALSE)
  }
  
  if (length(levels(group)) < 2) {
    stop("'group' must contain at least two groups.", call. = FALSE)
  }
  
  if (any(table(group) < 3)) {
    stop("Each group must contain at least three observations.", call. = FALSE)
  }
  
  if (!is.numeric(N) || N < 1) {
    stop("'N' must be a positive integer.", call. = FALSE)
  }
  
  N <- as.integer(N)
  
  if (!is.numeric(alpha) || alpha <= 0.5 || alpha > 1) {
    stop("'alpha' must be a numeric value in the interval (0.5, 1].",
         call. = FALSE)
  }
  
  if (!is.numeric(tau) || tau <= 0 || tau >= 1) {
    stop("'tau' must be a numeric value in the interval (0, 1).",
         call. = FALSE)
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  observed <- .RobHDMANOVA_stat(x = x, group = group, alpha = alpha,
                                tau = tau, cutoff = cutoff)
  
  TR_obs <- observed$TR
  perm_stats <- numeric(N)
  
  for (b in seq_len(N)) {
    perm_group <- sample(group)
    perm_stats[b] <- .RobHDMANOVA_stat(x = x, group = perm_group,
                                       alpha = alpha, tau = tau,
                                       cutoff = cutoff)$TR
    
    if (verbose && b %% 10 == 0) {
      message("Permutation ", b, " of ", N, " completed.")
    }
  }
  
  pval <- sum(perm_stats >= TR_obs) / N
  
  results <- list(
    Lambda = observed$Lambda,
    TR = TR_obs,
    p.value = pval,
    Permutations_TR = perm_stats,
    alpha = alpha,
    tau = tau,
    cutoff = cutoff,
    group.centers = observed$group.centers,
    weights = observed$weights,
    Test = "RobHDMANOVA"
  )
  
  class(results) <- c("MVTests", "list")
  return(results)
}


.RobHDMANOVA_stat <- function(x, group, alpha, tau, cutoff) {
  
  group <- as.factor(group)
  groups <- levels(group)
  g <- length(groups)
  n <- nrow(x)
  p <- ncol(x)
  
  group_centers_initial <- vector("list", g)
  group_sizes <- numeric(g)
  
  for (k in seq_len(g)) {
    xk <- x[group == groups[k], , drop = FALSE]
    group_sizes[k] <- nrow(xk)
    
    fit_k <- rrcov::CovMrcd(x = xk, alpha = alpha)
    group_centers_initial[[k]] <- as.numeric(fit_k@center)
  }
  
  centered_x <- x
  
  for (k in seq_len(g)) {
    idx <- which(group == groups[k])
    centered_x[idx, ] <- sweep(x[idx, , drop = FALSE], 2,
                               group_centers_initial[[k]], "-")
  }
  
  pooled_fit <- rrcov::CovMrcd(x = centered_x, alpha = alpha)
  WR <- as.matrix(pooled_fit@cov)
  delta <- as.numeric(pooled_fit@center)
  
  group_centers_preliminary <- vector("list", g)
  
  for (k in seq_len(g)) {
    group_centers_preliminary[[k]] <- group_centers_initial[[k]] + delta
  }
  
  inv_WR <- solve(WR)
  
  weights <- vector("list", g)
  group_centers <- vector("list", g)
  
  for (k in seq_len(g)) {
    xk <- x[group == groups[k], , drop = FALSE]
    
    diffs <- sweep(xk, 2, group_centers_preliminary[[k]], "-")
    rd <- sqrt(rowSums((diffs %*% inv_WR) * diffs))
    
    if (cutoff == "chisq") {
      crit <- sqrt(stats::qchisq(tau, df = p))
    } else {
      crit <- stats::median(rd) + stats::qnorm(tau) *
        stats::mad(rd, constant = 1.4826)
    }
    
    wk <- as.numeric(rd <= crit)
    
    if (sum(wk) == 0) {
      stop("All observations received zero weight in at least one group. ",
           "Try increasing 'tau' or 'alpha'.", call. = FALSE)
    }
    
    weights[[k]] <- wk
    group_centers[[k]] <- colSums(xk * wk) / sum(wk)
  }
  
  total_weight <- sum(unlist(weights))
  
  overall_center <- rep(0, p)
  
  for (k in seq_len(g)) {
    overall_center <- overall_center + sum(weights[[k]]) * group_centers[[k]]
  }
  
  overall_center <- overall_center / total_weight
  
  BR <- matrix(0, nrow = p, ncol = p)
  
  for (k in seq_len(g)) {
    diff_b <- matrix(group_centers[[k]] - overall_center, ncol = 1)
    BR <- BR + sum(weights[[k]]) * (diff_b %*% t(diff_b))
  }
  
  TR_mat <- WR + BR
  
  logdet_WR <- as.numeric(determinant(WR, logarithm = TRUE)$modulus)
  logdet_TR <- as.numeric(determinant(TR_mat, logarithm = TRUE)$modulus)
  
  Lambda <- exp(logdet_WR - logdet_TR)
  TR <- -log(Lambda)
  
  return(list(
    Lambda = Lambda,
    TR = TR,
    group.centers = group_centers,
    weights = weights
  ))
}