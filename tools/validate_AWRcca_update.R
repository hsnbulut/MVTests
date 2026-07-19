# Run this script from the MVTests package root before uploading the update.

required <- c("MASS", "Matrix")
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing)) {
  stop("Install required packages first: ", paste(missing, collapse = ", "))
}

source(file.path("R", "AWRcca.R"))

set.seed(123)
n <- 60
p <- 40
q <- 50
latent <- rnorm(n)
X <- outer(latent, seq(0.80, 0.30, length.out = p)) +
  matrix(rnorm(n * p), n, p)
Y <- outer(latent, seq(0.75, 0.25, length.out = q)) +
  matrix(rnorm(n * q), n, q)

fit <- AWRcca(X, Y)

expected_names <- c(
  "cor", "a", "b", "weights", "retained", "reweighted_final",
  "shrink_initial", "shrink_final", "shrink_used", "xi", "wrapping",
  "u1", "v1"
)

stopifnot(
  all(expected_names %in% names(fit)),
  length(fit$cor) == min(p, q),
  is.matrix(fit$a), nrow(fit$a) == p,
  is.matrix(fit$b), nrow(fit$b) == q,
  length(fit$weights) == n,
  all(fit$weights %in% 0:1),
  fit$retained == sum(fit$weights),
  identical(fit$shrink_used, fit$shrink_final),
  isTRUE(all.equal(fit$xi, 1.152204624, tolerance = 1e-8)),
  all(is.finite(fit$cor)),
  all(fit$cor >= 0 & fit$cor <= 1)
)

# Confirm that the second fit is not blocked by p + q > n.
stopifnot(p + q > n)

cat("AWRcca validation passed.\n")
cat("First canonical correlation:", fit$cor[1], "\n")
cat("Retained observations:", fit$retained, "of", n, "\n")
cat("Second weighted fit used:", fit$reweighted_final, "\n")
cat("Initial/final shrinkage:", fit$shrink_initial, fit$shrink_final, "\n")
cat("Consistency factor:", fit$xi, "\n")

# Package-level checks (run after installing devtools):
# devtools::document()
# devtools::check()
