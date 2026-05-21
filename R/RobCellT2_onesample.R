#' Cellwise Robust One-Sample Hotelling T^2 Test
#'
#' @description
#' Performs a cellwise robust one-sample Hotelling T^2 test based on the
#' cellwise minimum covariance determinant (cellMCD) estimator.
#'
#' @details
#' The function computes a robust Hotelling T^2 statistic by replacing the
#' classical sample mean vector and covariance matrix with the cellMCD location
#' and scatter estimates. The statistic is converted to an approximate F
#' statistic using the constants \code{d} and \code{q}. These constants can be
#' obtained by the \code{simRobCellT2_onesample()} function.
#'
#' @param data A numeric matrix or data frame.
#' @param mu0 The hypothesized mean vector under the null hypothesis.
#' @param d The scaling constant of the approximate F distribution.
#' @param q The second degree of freedom of the approximate F distribution.
#' @param alpha The cellMCD alpha parameter. Default is 0.75.
#' @param quant Quantile used in the cellMCD procedure. Default is 0.99.
#' @param crit Convergence criterion used in the cellMCD procedure. Default is 1e-04.
#' @param na.rm Logical. If TRUE, rows with missing values are removed. Default is TRUE.
#' @param ... Additional arguments passed to \code{cellWise::cellMCD()}.
#'
#' @return An object of class \code{MVTests} containing:
#' \item{T2}{The cellMCD-based robust Hotelling T^2 statistic.}
#' \item{Fval}{The approximate F statistic.}
#' \item{p.value}{The p-value based on the approximate F distribution.}
#' \item{mu}{The cellMCD location estimate.}
#' \item{S}{The cellMCD scatter estimate.}
#'
#' @references
#' Raymaekers, J. and Rousseeuw, P. J. (2024). The cellwise minimum covariance
#' determinant estimator. Journal of the American Statistical Association,
#' 119(548), 2610--2621.
#'
#' Willems, G., Pison, G., Rousseeuw, P. J., and Van Aelst, S. (2002).
#' A robust Hotelling test. Metrika, 55, 125--138.
#'
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#'
#' @examples
#' if (requireNamespace("MASS", quietly = TRUE) &&
#'     requireNamespace("cellWise", quietly = TRUE)) {
#'   set.seed(123)
#'   X <- MASS::mvrnorm(n = 50, mu = rep(0, 5), Sigma = diag(5))
#'
#'   const <- simRobCellT2_onesample(n = 50, p = 5, nrep = 50, seed = 123)
#'
#'   fit <- RobCellT2_onesample(
#'     data = X,
#'     mu0 = rep(0, 5),
#'     d = const$d,
#'     q = const$q
#'   )
#'
#'   fit$p.value
#' }
#'
#' @export
RobCellT2_onesample <- function(data, mu0, d, q,
                                alpha = 0.75, quant = 0.99, crit = 1e-04,
                                na.rm = TRUE, ...) {
  
  if (!requireNamespace("cellWise", quietly = TRUE)) {
    stop("The 'cellWise' package is required. Please install it first.",
         call. = FALSE)
  }
  
  data <- as.matrix(data)
  
  if (!is.numeric(data)) {
    stop("data must be a numeric matrix or data frame.", call. = FALSE)
  }
  
  if (na.rm) {
    data <- data[stats::complete.cases(data), , drop = FALSE]
  }
  
  n <- nrow(data)
  p <- ncol(data)
  
  if (n < 2) {
    stop("data must contain at least two complete observations.",
         call. = FALSE)
  }
  
  if (length(mu0) != p) {
    stop("The length of mu0 must be equal to the number of variables.",
         call. = FALSE)
  }
  
  if (!is.numeric(d) || length(d) != 1 || !is.finite(d) || d <= 0) {
    stop("d must be a positive finite numeric value.", call. = FALSE)
  }
  
  if (!is.numeric(q) || length(q) != 1 || !is.finite(q) || q <= 2) {
    stop("q must be a finite numeric value greater than 2.", call. = FALSE)
  }
  
  mu0 <- as.numeric(mu0)
  
  fit <- NULL
  invisible(utils::capture.output({
    fit <- tryCatch(
      cellWise::cellMCD(data, alpha = alpha, quant = quant, crit = crit, ...),
      error = function(e) NULL
    )
  }))
  
  if (is.null(fit)) {
    stop("cellMCD could not be computed for the supplied data.",
         call. = FALSE)
  }
  
  mu <- if (!is.null(fit$mu)) fit$mu else fit$center
  S  <- if (!is.null(fit$S))  fit$S  else fit$cov
  
  if (is.null(mu) || is.null(S)) {
    stop("The cellMCD output does not contain location and scatter estimates.",
         call. = FALSE)
  }
  
  mu <- as.numeric(mu)
  S  <- as.matrix(S)
  
  delta <- mu - mu0
  
  T2 <- tryCatch(
    as.numeric(n * t(delta) %*% solve(S) %*% delta),
    error = function(e) NA_real_
  )
  
  if (!is.finite(T2)) {
    stop("The robust Hotelling T^2 statistic could not be computed.",
         call. = FALSE)
  }
  
  Fval <- T2 / d
  pval <- stats::pf(Fval, df1 = p, df2 = q, lower.tail = FALSE)
  
  result <- list(
    Test = "RobCellT2_onesample",
    T2 = T2,
    statistic = T2,
    Fval = Fval,
    F = Fval,
    p.value = pval,
    pval = pval,
    df = c(p, q),
    parameter = c(n = n, p = p, d = d, q = q),
    method = "Cellwise robust one-sample Hotelling T^2 test",
    data.name = deparse(substitute(data)),
    mu0 = mu0,
    mu = mu,
    S = S,
    alpha.cellMCD = alpha,
    quant = quant,
    crit = crit
  )
  
  class(result) <- "MVTests"
  result
}