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
