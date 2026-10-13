# MVTests 2.3.6

* Added the X/mu0/Sigma/n/p/nper interface for one-sample cellMCD calibration.
  The null mean is required. When Sigma is omitted, it is estimated from raw X
  using the bias-corrected cellMCD scatter; a supplied Sigma takes precedence.
* Validate supplied dimensions against X and retain explicit n/p/Sigma mode.
  Return calibration dimensions, mean, covariance, moments and success/failure
  counts, with safeguards for finite-variance moment matching (q > 4).
* Preserve the former nrep/seed/data interface and fixed-center plug-in behavior
  in the deprecated simRobCellT2_onesample_legacy() compatibility function.
  Use set.seed() for reproducibility with the new interface.
* Updated examples, help and tests. The scaled-F reference remains an empirical
  approximation; the RobCellT2_onesample() statistic calculation is unchanged.

# MVTests 2.3.5

* Added named data-based null calibration using the cellMCD plug-in scatter,
  fixed null mean, matching complete-case handling, and recorded fit settings;
  retained the explicit-covariance positional API.

* Added Sigma-aware Monte Carlo calibration for the one-sample cellMCD test,
  with identity-model warnings, covariance validation, and explicit safeguards
  for finite-variance F moment matching.

# MVTests 2.3.4

* Updated `AWRcca()` so that score-space observation weights enter a second
  shrinkage-regularized canonical correlation fit. The final estimator no
  longer falls back to unregularized classical CCA merely because `p + q > n`.
* Replaced the Monte Carlo consistency correction with deterministic numerical
  integration.
* Added canonical directions and estimator diagnostics to the returned object
  while retaining the existing `cor`, `shrink_used`, `weights`, `u1`, and `v1`
  components for backward compatibility.
* Expanded input validation and numerical safeguards.
