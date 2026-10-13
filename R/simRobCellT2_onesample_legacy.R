#' Legacy Interface for One-Sample cellMCD Calibration
#'
#' @description
#' Compatibility interface for calls written before version 2.3.6. Prefer
#' \code{simRobCellT2_onesample(X, mu0, Sigma, n, p, nper)} for new code.
#' This interface retains the old n/p/nrep/seed/data arguments and the old
#' fixed-center plug-in design, and signals deprecation.
#' Computes the constants \code{d} and \code{q} required for the approximate
#' F distribution of the cellMCD-based robust one-sample Hotelling T^2 statistic.
#'
#' @param n The sample size in explicit-covariance mode; omit when supplying data.
#' @param p The number of variables in explicit-covariance mode; omit with data.
#' @param nrep The number of Monte Carlo replications. Default is 3000.
#' @param alpha The cellMCD alpha parameter. Default is 0.75.
#' @param quant Quantile used in the cellMCD procedure. Default is 0.99.
#' @param crit Convergence criterion used in the cellMCD procedure. Default is 1e-04.
#' @param seed Optional random seed.
#' @param Sigma Finite numeric positive-definite p-by-p covariance matrix under
#'   the null model. Defaults to \code{diag(p)}. Omitting it warns that the
#'   calibration assumes independence and may be inappropriate for correlated data.
#'   Supply either explicit n, p and Sigma or data, not both.
#' @param ... Additional arguments passed to \code{cellWise::cellMCD()}.
#' @param data Optional observed numeric matrix or data frame. Supplying data
#'   selects data-based plug-in calibration and derives n and p after row removal.
#' @param mu0 Finite numeric null mean vector of length p. Required with data;
#'   defaults to the zero vector in explicit-covariance mode.
#' @param na.rm Logical; remove incomplete rows when TRUE (default), matching
#'   \code{RobCellT2_onesample()}. With FALSE, missing cells are retained for
#'   cellMCD; infinite values are always rejected.
#'
#' @details
#' The scaled-F reference is an empirical, moment-matched approximation, not
#' a theoretically exact distribution. Calibration must match the test's sample
#' size, dimension, null covariance and cellMCD settings (including \code{...}).
#'
#' In data-based mode, cellMCD is fit to data centered at mu0 with
#' \code{fixedCenter = TRUE}, when supported by the installed API. Its
#' bias-corrected scatter \code{S} is the robust plug-in null covariance.
#' If fixedCenter is unavailable, a free-center fit of the centered data is used
#' with a warning and recorded in the settings. Failed fits never fall back to
#' classical covariance. The observed location estimate is not the simulation
#' mean: every sample is generated from \eqn{N_p(mu0, Sigma)} and refit to
#' calculate T^2 from the fitted location minus the fixed mu0. This is parametric
#' bootstrap/Monte Carlo calibration, not permutation calibration.
#'
#' All cellMCD settings agree between the plug-in and replicate fits except that
#' the replicate location is free. In data-based mode, do not set
#' \code{fixedCenter = TRUE} through \code{...}; that would fix the location
#' whose deviation the test measures. The plug-in covariance is an estimate,
#' not an exactly unbiased or guaranteed recovery of population covariance.
#' The approximation depends on that estimate and is not universally covariance
#' invariant. A supplied seed is set before the plug-in fit and simulations.
#'
#' Moment matching requires positive finite moments and a finite second degree
#' of freedom greater than 4. Incompatible or numerically unstable moments
#' produce an error; increasing \code{nrep} may help.
#'
#' @return A list with the following elements:
#' \item{d}{The scaling constant of the approximate F distribution.}
#' \item{q}{The second degree of freedom of the approximate F distribution.}
#' \item{mean.T2}{The Monte Carlo mean of the simulated T^2 statistics.}
#' \item{var.T2}{The Monte Carlo variance of the simulated T^2 statistics.}
#' \item{Sigma}{The covariance matrix used to generate calibration samples.}
#' \item{mode}{Either \code{explicit} or \code{data}.}
#' \item{n, p, mu0}{Effective sample size, dimension and simulation null mean.}
#' \item{n.omitted, na.rm}{Number of omitted rows and missing-data setting.}
#' \item{cellMCD.settings, covariance.settings}{Actual arguments, including API
#' defaults, for replicate and plug-in fits (the latter is NULL in explicit mode).}
#' \item{covariance.source, cellWise.version, seed}{Covariance source, installed
#' API version and supplied seed.}
#' \item{n.failed}{Number of failed or invalid replications.}
#' \item{n.success}{The number of successful Monte Carlo replications.}
#'
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#'
#' @examples
#' if (requireNamespace("MASS", quietly = TRUE) &&
#'     requireNamespace("cellWise", quietly = TRUE)) {
#'   mu0 <- c(1, -1)
#'   Sigma <- matrix(c(2, 0.4, 0.4, 1), 2)
#'   set.seed(123)
#'   X <- MASS::mvrnorm(30, mu0, Sigma)
#'   suppressWarnings(simRobCellT2_onesample_legacy(
#'     data = X, mu0 = mu0, nrep = 50, seed = 1))
#'   suppressWarnings(simRobCellT2_onesample_legacy(
#'     30, 2, nrep = 50, seed = 1, Sigma = Sigma))
#' }
#'
#' @export
simRobCellT2_onesample_legacy <- function(n, p, nrep = 3000,
                                          alpha = 0.75, quant = 0.99, crit = 1e-04,
                                          seed = NULL, Sigma = diag(p), ...,
                                          data = NULL, mu0 = NULL, na.rm = TRUE) {
  .Deprecated("simRobCellT2_onesample", package = "MVTests",
              msg = "The legacy calibration API is deprecated. Migrate to simRobCellT2_onesample(X, mu0, Sigma, n, p, nper); set.seed() replaces seed. The new plug-in covariance uses raw X, while this compatibility interface retains the old fixed-center design.")
  call <- match.call()
  call[[1L]] <- .simRobCellT2_onesample_legacy_impl
  eval(call, envir = parent.frame())
}

.simRobCellT2_onesample_legacy_impl <- function(n, p, nrep = 3000,
                                   alpha = 0.75, quant = 0.99, crit = 1e-04,
                                   seed = NULL, Sigma = diag(p), ...,
                                   data = NULL, mu0 = NULL, na.rm = TRUE) {

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required. Please install it first.",
         call. = FALSE)
  }

  if (!requireNamespace("cellWise", quietly = TRUE)) {
    stop("The 'cellWise' package is required. Please install it first.",
         call. = FALSE)
  }

  valid.integer <- function(x, minimum) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= minimum && x <= .Machine$integer.max && x == floor(x)
  }
  if (!valid.integer(nrep, 1)) {
    stop("nrep must be a positive finite integer.", call. = FALSE)
  }
  nrep <- as.integer(nrep)
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("na.rm must be TRUE or FALSE.", call. = FALSE)
  }
  Sigma.omitted <- missing(Sigma)
  mode <- if (is.null(data)) "explicit" else "data"
  n.omitted <- 0L
  if (mode == "data") {
    if (!missing(n) || !missing(p) || !missing(Sigma)) {
      stop("Supply data and mu0 OR explicit n, p and Sigma, not both modes.",
           call. = FALSE)
    }
    data <- tryCatch(as.matrix(data), error = function(e) NULL)
    if (is.null(data) || !is.numeric(data) || length(dim(data)) != 2L) {
      stop("data must be a numeric matrix or data frame.", call. = FALSE)
    }
    if (na.rm) {
      keep <- stats::complete.cases(data)
      n.omitted <- sum(!keep)
      data <- data[keep, , drop = FALSE]
    }
    if (any(!is.finite(data) & !is.na(data))) {
      stop("data must contain finite values in all non-missing cells.",
           call. = FALSE)
    }
    n <- nrow(data)
    p <- ncol(data)
    if (is.null(mu0)) {
      stop("mu0 is required for data-based calibration.", call. = FALSE)
    }
  } else if (missing(n) || missing(p)) {
    stop("Supply n and p for explicit-covariance calibration, or data and mu0.",
         call. = FALSE)
  }
  if (!valid.integer(n, 2)) {
    stop("n must be a finite integer greater than 1 (after complete-case handling).",
         call. = FALSE)
  }
  if (!valid.integer(p, 2)) {
    stop("p must be a finite integer greater than 1.", call. = FALSE)
  }
  n <- as.integer(n)
  p <- as.integer(p)
  if (is.null(mu0)) mu0 <- rep(0, p)
  if (!is.numeric(mu0) || !is.null(dim(mu0)) ||
      length(mu0) != p || any(!is.finite(mu0))) {
    stop("mu0 must be a finite numeric vector of length p.", call. = FALSE)
  }
  mu <- as.numeric(mu0)

  validate.Sigma <- function(Sigma, label = "Sigma") {
    if (!is.matrix(Sigma) || !is.numeric(Sigma) ||
        !identical(dim(Sigma), c(p, p)) || any(!is.finite(Sigma))) {
      stop(paste(label, "must be a finite numeric p-by-p matrix."), call. = FALSE)
    }
    if (!isSymmetric(Sigma, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE)) {
      stop(paste(label, "must be symmetric within numerical tolerance."),
           call. = FALSE)
    }
    Sigma <- Sigma / 2 + t(Sigma) / 2
    if (is.null(tryCatch(chol(Sigma), error = function(e) NULL))) {
      stop(paste(label, "must be positive definite."), call. = FALSE)
    }
    Sigma
  }

  # Resolve partial/positional cellMCD arguments once, catching duplicates early.
  fit.args <- tryCatch(
    as.list(match.call(cellWise::cellMCD,
      as.call(c(list(quote(cellMCD), X = quote(.X)),
                list(alpha = alpha, quant = quant, crit = crit), list(...)))))[-1L],
    error = function(e) stop(paste("Invalid cellMCD arguments:", conditionMessage(e)),
                             call. = FALSE)
  )
  fit.args$X <- NULL
  api <- formals(cellWise::cellMCD)
  for (name in setdiff(names(api), c("X", "...", names(fit.args)))) {
    fit.args[name] <- list(eval(api[[name]], envir = environment(cellWise::cellMCD)))
  }
  covariance.settings <- NULL
  if (mode == "data") {
    if ("fixedCenter" %in% names(fit.args) &&
        !identical(fit.args$fixedCenter, FALSE)) {
      stop("Data-based calibration requires free replicate locations: fixedCenter must be FALSE in ... .",
           call. = FALSE)
    }
    covariance.settings <- fit.args
    if ("fixedCenter" %in% names(api)) {
      covariance.settings$fixedCenter <- TRUE
    } else {
      warning("Installed cellWise API does not support fixedCenter; the centered-data plug-in fit uses a free center.",
              call. = FALSE)
    }
  }
  if (!is.null(seed)) set.seed(seed)
  if (mode == "data") {
    fit0 <- NULL
    invisible(utils::capture.output({
      fit0 <- tryCatch(
        do.call(cellWise::cellMCD,
                c(list(X = sweep(data, 2L, mu, "-")), covariance.settings)),
        error = function(e) stop(paste("cellMCD plug-in null covariance estimation failed:",
                                      conditionMessage(e)), call. = FALSE)
      )
    }))
    Sigma <- validate.Sigma(fit0$S, "cellMCD plug-in covariance S")
  } else {
    if (missing(Sigma)) {
      warning("Sigma was omitted: d and q are calibrated under the independence/identity-covariance model and may not be appropriate for correlated data.",
              call. = FALSE)
    }
    Sigma <- validate.Sigma(Sigma)
  }

  T2_sim <- rep(NA_real_, nrep)
  failure.messages <- rep("", nrep)

  for (i in seq_len(nrep)) {

    X0 <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)

    fit <- NULL
    invisible(utils::capture.output({
      fit <- tryCatch(
        do.call(cellWise::cellMCD, c(list(X = X0), fit.args)),
        error = function(e) {
          failure.messages[i] <<- paste("cellMCD fit failed:", conditionMessage(e))
          NULL
        }
      )
    }))

    if (is.null(fit)) {
      if (!nzchar(failure.messages[i])) failure.messages[i] <- "cellMCD returned no fit."
      next
    }

    mu_hat <- if (!is.null(fit$mu)) fit$mu else fit$center
    S_hat  <- if (!is.null(fit$S))  fit$S  else fit$cov

    if (is.null(mu_hat) || is.null(S_hat)) {
      failure.messages[i] <- "cellMCD output lacks location or scatter."
      next
    }

    if (!is.numeric(mu_hat) || length(mu_hat) != p || any(!is.finite(mu_hat))) {
      failure.messages[i] <- "cellMCD location is invalid."
      next
    }
    delta <- as.numeric(mu_hat) - mu
    S_hat <- tryCatch(validate.Sigma(as.matrix(S_hat)), error = function(e) {
      failure.messages[i] <<- paste("cellMCD scatter is invalid:", conditionMessage(e))
      NULL
    })
    if (is.null(S_hat)) next

    T2_sim[i] <- tryCatch(
      as.numeric(n * t(delta) %*% solve(S_hat) %*% delta),
      error = function(e) {
        failure.messages[i] <<- paste("T2 calculation failed:", conditionMessage(e))
        NA_real_
      }
    )
    if ((!is.finite(T2_sim[i]) || T2_sim[i] <= 0) && !nzchar(failure.messages[i])) {
      failure.messages[i] <- "T2 must be positive and finite."
    }
  }

  failed <- !is.finite(T2_sim) | T2_sim <= 0
  failure.messages <- stats::setNames(failure.messages[failed], which(failed))
  T2_sim <- T2_sim[!failed]

  n.success <- length(T2_sim)

  if (n.success < 10) {
    stop(sprintf("Too few successful cellMCD replications (%d/%d). Increase nrep or check the data dimension.",
                 n.success, nrep),
         call. = FALSE)
  }

  muT2 <- mean(T2_sim)
  varT2 <- stats::var(T2_sim)

  if (!is.finite(muT2) || muT2 <= 0 ||
      !is.finite(varT2) || varT2 <= 0) {
    stop("Monte Carlo moments must be positive and finite. Increase nrep.",
         call. = FALSE)
  }
  moment.ratio <- (varT2 / muT2^2) * (p / 2)
  if (!is.finite(moment.ratio) ||
      moment.ratio - 1 <= sqrt(.Machine$double.eps) * max(1, moment.ratio)) {
    stop("Monte Carlo moments cannot stably match a finite-variance F distribution (q > 4). Increase nrep or check the calibration model.",
         call. = FALSE)
  }
  q <- (p + 2) / (moment.ratio - 1) + 4
  d <- muT2 * (q - 2) / q

  if (!is.finite(d) || !is.finite(q) || d <= 0 || q <= 4) {
    stop("Moment matching produced invalid d or q; a finite-variance F distribution requires d > 0 and q > 4. Increase nrep.",
         call. = FALSE)
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
    crit = crit,
    Sigma = Sigma,
    mode = mode,
    n = n,
    p = p,
    mu0 = mu,
    na.rm = na.rm,
    n.omitted = n.omitted,
    n.failed = nrep - n.success,
    failure.messages = failure.messages,
    seed = seed,
    cellWise.version = as.character(utils::packageVersion("cellWise")),
    covariance.source = if (mode == "data") "cellMCD S" else if (Sigma.omitted)
      "default identity" else "supplied Sigma",
    cellMCD.settings = fit.args,
    covariance.settings = covariance.settings
  )
}
