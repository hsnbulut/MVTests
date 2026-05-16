% Generated manually for MVTests
\name{CellMCDT2}
\alias{CellMCDT2}
\title{Cellwise Robust Permutation Hotelling T^2 Test for Two Independent Samples}
\usage{
CellMCDT2(X1, X2, B = 100, alpha = 0.75, quant = 0.99,
  crit = 1e-04, noCits = 100, lmin = 1e-04,
  seed = NULL, correction = TRUE, verbose = FALSE)
}
\arguments{
\item{X1}{A numeric matrix or data frame containing the observations from the first group.}
\item{X2}{A numeric matrix or data frame containing the observations from the second group.}
\item{B}{The number of random permutations used to approximate the null distribution.}
\item{alpha}{The cellMCD trimming parameter. In each column, at least \code{n * alpha} cells must remain unflagged.}
\item{quant}{Cutoff probability used by \code{cellWise::cellMCD} to flag outlying cells.}
\item{crit}{Convergence tolerance used by \code{cellWise::cellMCD}.}
\item{noCits}{Maximum number of C-steps used by \code{cellWise::cellMCD}.}
\item{lmin}{Lower bound on the eigenvalues of the estimated covariance matrix on the standardized data.}
\item{seed}{Optional integer used to set the random seed for the permutation procedure.}
\item{correction}{Logical. If \code{TRUE}, the permutation p-value is calculated as \code{(sum(T_perm >= T_obs) + 1) / (B + 1)}.}
\item{verbose}{Logical. If \code{TRUE}, progress information is printed during the permutation procedure.}
}
\description{
Performs the cellMCD-based robust permutation Hotelling T^2 test for comparing the mean vectors of two independent multivariate samples in the presence of cellwise outliers.
}
\details{
For each group, the cellMCD location vector and covariance matrix are computed. These estimates are then used to construct a robust pooled covariance matrix and the observed robust Hotelling T^2 statistic. Since the finite-sample null distribution of the statistic is unknown, the p-value is obtained by randomly permuting the group labels.
}
\value{
A list of class \code{MVTests} with elements including \code{T2}, \code{p.value}, \code{Permutations_T2}, \code{B}, \code{alpha}, \code{quant}, sample sizes, group centers, group covariance matrices, and the robust pooled covariance matrix.
}
\references{
Raymaekers, J. and Rousseeuw, P. J. (2024). The cellwise minimum covariance determinant estimator. \emph{Journal of the American Statistical Association}, 119(548), 2610--2621.

Bulut, H. and Esmeray, M. A cellwise robust Hotelling test for two-sample comparisons (Unpublished).
}
\author{Hasan BULUT \email{hasan.bulut@omu.edu.tr}}
\examples{
if (requireNamespace("mvtnorm", quietly = TRUE) &&
    requireNamespace("cellWise", quietly = TRUE)) {
  set.seed(123)
  x1 <- mvtnorm::rmvnorm(n = 30, mean = rep(0, 5), sigma = diag(5))
  x2 <- mvtnorm::rmvnorm(n = 30, mean = rep(0, 5), sigma = diag(5))
  fit <- CellMCDT2(X1 = x1, X2 = x2, B = 19, seed = 123)
  fit$p.value
}
}
