#' Monte Carlo Simulation for d and q Constants of RobCellT2_onesample
#'
#' @description
#' Computes the constants \code{d} and \code{q} required for the approximate
#' F distribution of the cellMCD-based robust one-sample Hotelling T^2 statistic.
#'
#' @param n The sample size.
#' @param p The number of variables.
#' @param nrep The number of Monte Carlo replications. Default is 3000.
#' @param alpha The cellMCD alpha parameter. Default is 0.75.
#' @param quant Quantile used in the cellMCD procedure. Default is 0.99.
#' @param crit Convergence criterion used in the cellMCD procedure. Default is 1e-04.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to \code{cellWise::cellMCD()}.
#'
#' @return A list with the following elements:
#' \item{d}{The scaling constant of the approximate F distribution.}
#' \item{q}{The second degree of freedom of the approximate F distribution.}
#' \item{mean.T2}{The Monte Carlo mean of the simulated T^2 statistics.}
#' \item{var.T2}{The Monte Carlo variance of the simulated T^2 statistics.}
#' \item{n.success}{The number of successful Monte Carlo replications.}
#'
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#'
#' @examples
#' if (requireNamespace("MASS", quietly = TRUE) &&
#'     requireNamespace("cellWise", quietly = TRUE)) {
#'   simRobCellT2_onesample(n = 50, p = 5, nrep = 50, seed = 123)
#' }
#'
#' @export
simRobCellT2_onesample <- function(n, p, nrep = 3000,
                                   alpha = 0.75, quant = 0.99, crit = 1e-04,
                                   seed = NULL, ...) {
  
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required. Please install it first.",
         call. = FALSE)
  }
  
  if (!requireNamespace("cellWise", quietly = TRUE)) {
    stop("The 'cellWise' package is required. Please install it first.",
         call. = FALSE)
  }
  
  if (!is.numeric(n) || length(n) != 1 || n <= 1) {
    stop("n must be a numeric value greater than 1.", call. = FALSE)
  }
  
  if (!is.numeric(p) || length(p) != 1 || p <= 1) {
    stop("p must be a numeric value greater than 1.", call. = FALSE)
  }
  
  if (!is.numeric(nrep) || length(nrep) != 1 || nrep < 1) {
    stop("nrep must be a positive integer.", call. = FALSE)
  }
  
  n <- as.integer(n)
  p <- as.integer(p)
  nrep <- as.integer(nrep)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  mu <- rep(0, p)
  sigma <- diag(p)
  
  T2_sim <- rep(NA_real_, nrep)
  
  for (i in seq_len(nrep)) {
    
    X0 <- MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
    
    fit <- NULL
    invisible(utils::capture.output({
      fit <- tryCatch(
        cellWise::cellMCD(X0, alpha = alpha, quant = quant, crit = crit, ...),
        error = function(e) NULL
      )
    }))
    
    if (is.null(fit)) {
      next
    }
    
    mu_hat <- if (!is.null(fit$mu)) fit$mu else fit$center
    S_hat  <- if (!is.null(fit$S))  fit$S  else fit$cov
    
    if (is.null(mu_hat) || is.null(S_hat)) {
      next
    }
    
    delta <- as.numeric(mu_hat) - mu
    S_hat <- as.matrix(S_hat)
    
    T2_sim[i] <- tryCatch(
      as.numeric(n * t(delta) %*% solve(S_hat) %*% delta),
      error = function(e) NA_real_
    )
  }
  
  T2_sim <- T2_sim[is.finite(T2_sim) & T2_sim > 0]
  
  n.success <- length(T2_sim)
  
  if (n.success < 10) {
    stop("Too few successful cellMCD replications. Increase nrep or check the data dimension.",
         call. = FALSE)
  }
  
  muT2 <- mean(T2_sim)
  varT2 <- stats::var(T2_sim)
  
  q <- (1 / (((varT2 / muT2^2) * (p / 2)) - 1)) * (p + 2) + 4
  d <- muT2 * (q - 2) / q
  
  if (!is.finite(d) || !is.finite(q) || d <= 0 || q <= 2) {
    warning("The estimated d or q value is not stable. Consider increasing nrep.")
  }
  
  list(
    d = d,
    q = q,
    mean.T2 = muT2,
    var.T2 = varT2,
    n.success = n.success,
    nrep = nrep,
    alpha.cellMCD = alpha,
    quant = quant,
    crit = crit
  )
}