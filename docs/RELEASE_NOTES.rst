MotifSurvey 0.2.1
=================

This release adds a reproducible sampling-design case study for network closure.
Three independent-vertex allocations use department-roster information and the
same expected vertex budget. Separate pilot and evaluation samples distinguish
design selection from evaluation; every candidate allocation is reported.

The example computes the linear contrast 3 * triangles - 0.25 * wedges with
its estimated variance, retaining the wedge/triangle covariance. It does not
provide confidence intervals or validate the original collection mechanism.

Changes
-------
* Add the example, frozen protocol, exact reference-population covariance
  oracle, full realization records and reproduction instructions.
* Add five tests, including all four-node graphs checked against literal
  motif pairs and an exhaustive sampling-distribution check.
* Centralize version reporting so package, API and CLI versions agree.
* Make failed numerical stress checks return a failing process status.
* Add a portable installed-distribution checker and a release workflow that
  requires the Windows/Linux/macOS matrix and minimum-dependency check.

The estimation algorithm is unchanged from 0.2.0. Its core-file differences
only import and report the shared version. Historical validation records keep
their original version labels. CI artifacts identify their actual platform,
interpreter, dependency versions, commit and run.

The software supports a fixed simple undirected graph, independently selected
vertices with known positive probabilities, and complete induced edges. It
still enumerates every edge, wedge and triangle and enforces resource budgets.
External adoption and independent human statistical review are not claimed.
