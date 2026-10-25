API and input format
====================

``estimate(nodes, edges, probabilities, *, design, complete_induced,
max_instances=200000, max_subsets=1000000)`` returns a JSON-compatible dict.
Callers must pass ``design="independent_node"`` and ``complete_induced=True``.

``nodes`` contains distinct nonempty string IDs, including selected isolates.
``edges`` contains unique unordered pairs of these IDs, without self-loops.
``probabilities`` maps exactly those IDs to finite numeric values in (0, 1].
The complete population and its unselected vertices are not supplied.
Positive inclusion probability for every population vertex is an external
design requirement. Numeric strings and Boolean probabilities are rejected.

The observed graph must be simple, undirected and completely induced by the
selected vertices. Any conversion from a directed graph changes the estimand
and must be described by the researcher before calling the package.

Output fields
-------------
* ``software``, ``version``, ``design`` and ``estimand`` identify the analysis.
* ``motif_order`` is always edges, wedges, triangles. All vector and matrix
  entries follow this order. Wedges are centered unordered neighbor pairs.
* ``observed_counts`` contains unweighted counts in the observed graph.
* ``estimated_totals`` contains the three Horvitz-Thompson estimates.
* ``estimated_covariance`` is the raw symmetric 3 by 3 design-covariance
  estimate. It may be indefinite on a particular sample.
* ``warnings`` flags zero counts, indefinite covariance, negative diagonal
  from cancellation, and the absence of a certified normal approximation.
* ``diagnostics`` records graph size, enumerated units, subset keys, the
  smallest covariance eigenvalue and absolute covariance-term sums.
* ``input_sha256`` hashes sorted normalized inputs and design declarations.
  It identifies the input, not an external verification of the sampling design.

The absolute-term sum can help identify cancellation; it is not a formal
error bound. Hash stability is specified for the normalized JSON produced
by this version. No confidence interval or covariance correction is applied.

Resources and errors
--------------------
Both resource limits are positive integers. Count checks stop before motif
enumeration when the edge-plus-wedge lower bound already exceeds the budget.
Triangle and subset counts are checked during accumulation. Raising limits
is explicit; there is no automatic fallback to approximate counting.

Invalid declarations, malformed values and nonfinite computations raise
``ValueError`` (some structurally invalid Python objects raise ``TypeError``).
``ResourceLimitError`` is a ``ValueError`` subclass. The CLI catches these
and file errors, reports the problem to stderr and exits with code 2.
Successful output exits with code 0. Existing output paths are never replaced.

Development
-----------
Keep normalization, counting, weighting and covariance behavior covered by
the existing tests. Changes to the statistical core require repeating the
direct-pair and exact-expectation checks, recording the new core hash and
updating benchmark version labels. Do not reuse timings from a different
version as new measurements. Source and wheel installations must both pass.
