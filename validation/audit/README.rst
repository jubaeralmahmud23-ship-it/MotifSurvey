Mathematics, code and result audit: 0.2.2
========================================

This is an internal computational audit, not independent human peer review.
The exact-fraction implementation enumerates triples and centers independently
of the package and evaluates literal ordered motif pairs. For every four-node
graph, two probability regimes and all 16 vertex subsets, it checks the package
totals and covariance. Exact expectations equal the population totals and true
design covariance in all 128 sampling laws (2,048 sample calculations).

All 2,400 archived pilot/evaluation contrasts and variances are recomputed from
the saved three totals and covariance matrices. MSE, variance summaries, pilot
selection, negative variance counts and indefinite-matrix counts agree. These
are algebraic rechecks of archived records, not new independent samples.
The largest discrepancy between mean estimated and exact variance is 3.035
Monte Carlo standard errors. The finite study therefore does not establish
tight variance calibration or interval coverage. The unbiasedness result
rests on the derivation and exact checks, not a claim that every simulation
mean is close to its expectation.

Run from an extracted source distribution after installing it::

    python validation/audit/run_audit.py --output new-audit-results.json

The output must be a new path. results.json records the executed local run;
its teaching_examples field contains the worked edge and triangle examples.
The platform workflow reruns the audit and retains each environment's report.

The code audit found that NumPy Boolean/complex scalars and byte strings could
be converted to probabilities, and huge integers could raise OverflowError.
Version 0.2.2 rejects these inputs with ValueError. Threshold validation in
the closure example was tightened similarly. Three regression tests bring
the core suite to 31 tests; the companion has four tests. The estimation body
from adjacency construction onward is byte-identical, after normalizing line
endings, to version 0.2.1. The experiments used valid real probabilities and
their reported values are unaffected. Historical timings retain their versions.

Reference and claim checks are recorded in docs/REFERENCE_AUDIT.rst.
