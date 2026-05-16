#' Cellwise Robust Permutation Hotelling T^2 Test for Two Independent Samples
#'
#' Performs a cellMCD-based robust two-sample Hotelling T^2 test for comparing
#' the mean vectors of two independent multivariate samples. The p-value is
#' obtained by a permutation procedure.
#'
#' @param X1 A numeric matrix or data frame for the first group.
#' @param X2 A numeric matrix or data frame for the second group.
#' @param B Number of permutations. Default is 999.
#' @param alpha The cellMCD alpha parameter. Default is 0.75.
#' @param quant Quantile used in the cellMCD procedure. Default is 0.99.
#' @param crit Convergence criterion used in the cellMCD procedure. Default is 1e-04.
#' @param seed Optional random seed.
#' @param na.rm Logical. If TRUE, rows with missing values are removed. Default is TRUE.
#' @param ... Additional arguments passed to \code{cellWise::cellMCD()}.
#'
#' @return An object of class \code{MVTests} containing the test statistic,
#' permutation p-value, number of successful permutations, and related information.
#'
#' @references
#' Raymaekers, J. and Rousseeuw, P. J. (2024). The cellwise minimum covariance
#' determinant estimator. Journal of the American Statistical Association,
#' 119(548), 2610--2621.
#'
#' Bulut, H. and Esmeray, M. A cellwise robust Hotelling test for two-sample comparisons (Unpublished).
#'
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#' 
#' @examples
#' if (requireNamespace("mvtnorm", quietly = TRUE) &&
#'     requireNamespace("cellWise", quietly = TRUE)) {
#'   set.seed(123)
#'   x1 <- mvtnorm::rmvnorm(n = 30, mean = rep(0, 5), sigma = diag(5))
#'   x2 <- mvtnorm::rmvnorm(n = 30, mean = rep(0, 5), sigma = diag(5))
#'
#'   fit <- CellMCDT2(X1 = x1, X2 = x2, B = 9, seed = 123)
#'   fit$p.value
#' }
#'
#' @export


CellMCDT2 <- function(X1, X2, B = 999, alpha = 0.75, quant = 0.99,
                      crit = 1e-04, seed = NULL, na.rm = TRUE, ...) {
  
  if (!requireNamespace("cellWise", quietly = TRUE)) {
    stop("The 'cellWise' package is required. Please install it first.",
         call. = FALSE)
  }
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  if (!is.numeric(X1) || !is.numeric(X2)) {
    stop("X1 and X2 must be numeric matrices or data frames.", call. = FALSE)
  }
  
  if (ncol(X1) != ncol(X2)) {
    stop("X1 and X2 must have the same number of variables.", call. = FALSE)
  }
  
  if (na.rm) {
    X1 <- X1[stats::complete.cases(X1), , drop = FALSE]
    X2 <- X2[stats::complete.cases(X2), , drop = FALSE]
  }
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p  <- ncol(X1)
  
  if (n1 < 2 || n2 < 2) {
    stop("Each group must contain at least two observations.", call. = FALSE)
  }
  
  if (B < 1) {
    stop("B must be a positive integer.", call. = FALSE)
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  .quiet_cellMCD <- function(X) {
    tmp <- NULL
    invisible(utils::capture.output({
      tmp <- tryCatch(
        cellWise::cellMCD(X, alpha = alpha, quant = quant, crit = crit, ...),
        error = function(e) NULL
      )
    }))
    tmp
  }
  
  .extract_cellmcd <- function(fit) {
    if (is.null(fit)) {
      return(NULL)
    }
    
    mu <- NULL
    S  <- NULL
    
    if (!is.null(fit$mu)) {
      mu <- fit$mu
    } else if (!is.null(fit$center)) {
      mu <- fit$center
    }
    
    if (!is.null(fit$S)) {
      S <- fit$S
    } else if (!is.null(fit$cov)) {
      S <- fit$cov
    }
    
    if (is.null(mu) || is.null(S)) {
      return(NULL)
    }
    
    list(mu = as.numeric(mu), S = as.matrix(S))
  }
  
  .stat_fun <- function(A, Bmat) {
    nA <- nrow(A)
    nB <- nrow(Bmat)
    
    fitA <- .extract_cellmcd(.quiet_cellMCD(A))
    fitB <- .extract_cellmcd(.quiet_cellMCD(Bmat))
    
    if (is.null(fitA) || is.null(fitB)) {
      return(NA_real_)
    }
    
    Cpool <- ((nA - 1) * fitA$S + (nB - 1) * fitB$S) / (nA + nB - 2)
    
    diff <- fitA$mu - fitB$mu
    
    stat <- tryCatch(
      as.numeric((nA * nB / (nA + nB)) *
                   t(diff) %*% solve(Cpool) %*% diff),
      error = function(e) NA_real_
    )
    
    stat
  }
  
  T.obs <- .stat_fun(X1, X2)
  
  if (is.na(T.obs)) {
    stop("The observed cellMCD-based statistic could not be computed.",
         call. = FALSE)
  }
  
  Xpool <- rbind(X1, X2)
  n <- nrow(Xpool)
  
  T.perm <- rep(NA_real_, B)
  
  for (b in seq_len(B)) {
    ind <- sample.int(n, size = n1, replace = FALSE)
    X1b <- Xpool[ind, , drop = FALSE]
    X2b <- Xpool[-ind, , drop = FALSE]
    
    T.perm[b] <- .stat_fun(X1b, X2b)
  }
  
  T.perm.valid <- T.perm[is.finite(T.perm)]
  B.success <- length(T.perm.valid)
  
  if (B.success == 0) {
    stop("No successful permutation statistics could be computed.",
         call. = FALSE)
  }
  
  p.value <- (sum(T.perm.valid >= T.obs) + 1) / (B.success + 1)
  
  result <- list(
    statistic = T.obs,
    p.value = p.value,
    method = "Cellwise robust permutation Hotelling T^2 test",
    data.name = paste(deparse(substitute(X1)), "and", deparse(substitute(X2))),
    parameter = c(
      n1 = n1,
      n2 = n2,
      p = p,
      B = B,
      successful.permutations = B.success
    ),
    permutation.statistics = T.perm,
    alpha.cellMCD = alpha,
    quant = quant,
    crit = crit
  )
  
  class(result) <- "MVTests"
  result
}