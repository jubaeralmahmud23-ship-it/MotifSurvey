Change log
==========

0.2.2 - input validation audit
------------------------------
* Reject NumPy Boolean and complex probabilities, text/byte values and
  non-scalar objects before float conversion. Keep finite real scalars,
  Decimal values and Fraction values supported.
* Convert overflowing real probabilities to a documented ValueError; the
  CLI reports exit code 2 without a traceback.
* Apply corresponding threshold checks in the closure-profile example.
* Add three probability regression tests (31 package tests in total).
* Preserve counting, weighting and covariance arithmetic for valid inputs.
  Earlier benchmark and study records retain their original version labels.

Research companion - 25 October 2026
------------------------------------
* Add a reusable closure-profile example, explicit tool-selection guidance,
  all 27 budget/allocation/threshold sensitivity settings and retrospective
  reanalysis of the original evaluation samples.
* Execute the native R survey adapter on the same three sampled research
  graphs and compare nine contrast/variance pairs.
* Add four companion checks, kept separate from the 28 release tests.
* Preserve the v0.2.1 tag, distributions and installed estimator unchanged.

0.2.1 - research case and release checks
--------------------------------------
* Add a department-roster sampling-design example for the linear closure
  contrast, with a frozen protocol, separate pilot/evaluation simulations,
  exact reference-population covariance and auditable raw results.
* Add five tests for the new mathematical example and version consistency.
* Centralize package/result version reporting; estimation arithmetic unchanged.
* Fail the stress-check process if any numerical or rejection case fails.
* Exercise installed wheels and source distributions in a portable checker;
  configure a 12-job platform/interpreter matrix and minimum-dependency job.
* Gate tagged publication on successful checks. Configuration alone is not
  evidence that a hosted job or public release has occurred.

0.2.0 - local release prepared 25 October 2026
--------------------------------------------
* Reduce subset covariance contributions in NumPy blocks of 512 keys.
* Preserve raw covariance and expose numerical and inference diagnostics.
* Add native R comparisons, high-precision stress cases, sparse baselines,
  controlled public-network thinning and a large-star closed-form check.
* Package source, wheel, examples, derivation and reproducible evidence.
* Select NetworkX 3.4.2 for the validation extra on Python 3.10 and 3.6.1
  on newer Python versions. This fixes the original validation-extra pin.

The statistical core is byte-identical to the validated 0.2.0 prototype.
Release preparation changes documentation and packaging, not the estimators.
This version has not yet been uploaded or tagged in a public repository.

0.1.0 - historical local prototype
----------------------------------
* Initial totals, subset covariance, validation and command-line interface.
* The frozen core and original measurements remain in validation/followup
  to make the version comparison auditable.
