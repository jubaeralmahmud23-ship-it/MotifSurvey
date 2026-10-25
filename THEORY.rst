Statistical specification and audit
==================================

This is an internal mathematical audit, not independent human peer review.
The estimator and covariance are classical sampling theory. The 2024 revision
of Bhattacharya, Das and Mukherjee gives the uniform-probability overlap
variance explicitly (Proposition 2.4), before HT rescaling by p^(-2k)
for a k-vertex motif: https://arxiv.org/abs/2011.03026v2.

Population and observations
---------------------------
The population is a fixed finite simple undirected graph. Each vertex v is
selected independently with known probability p(v) in (0,1]. The observed
graph contains all edges between selected vertices. Unobserved vertices need
not be listed in the input, but positive inclusion probability for every
population vertex is an external design assumption. Estimated probabilities,
fixed-size sampling, ego samples and missing observed edges are not supported.

A motif unit is a particular edge, centered unordered neighbor pair (wedge),
or triangle. Thus three wedge units share the same vertices in a triangle;
they remain distinct population units. Let pi(M) be the product of the vertex
probabilities in unit M, and I(M) its observation indicator. The count estimate
for type t is sum over observed M of type t of 1/pi(M).

Unbiased totals and covariance
------------------------------
E[I(M)] = pi(M), so each population motif contributes one in expectation.
For two motif units M and N, pi(M,N) is the product over their union.
The observed covariance contribution is

  I(M) I(N) / [pi(M) pi(N)] * [1 - pi(M) pi(N)/pi(M,N)].

Its expectation is [pi(M,N) - pi(M) pi(N)]/[pi(M) pi(N)],
the corresponding summand of the true design covariance. Independent vertex
selection simplifies the bracket to 1 minus the product of probabilities
over the intersection. Disjoint motif units therefore contribute zero.
All pairs are ordered, including diagonal pairs and distinct wedge units
with identical vertex sets. Cross-type entries are computed at the same time.

Matrix-free computation
------------------------
For each nonempty vertex subset A contained in an observed motif, let g(A,t)
be the sum of 1/pi(M) over type-t motifs containing A. Inclusion-exclusion
expands one minus the intersection product into

  sum over nonempty A in the intersection:
      (-1)^(|A|+1) product over v in A of (1-p(v)).

Consequently the covariance matrix is a signed sum of g(A)g(A)^T with these
coefficients. With motifs of at most three vertices, each motif updates at
most seven subset keys. This avoids storage of a motif-by-motif matrix.
It still enumerates every motif instance; it does not make dense-graph motif
enumeration constant-time or remove degree-dependent growth in wedge counts.

Version 0.2.0 computes this sum in blocks of at most 512 keys. The block
matrix multiplication is algebraically identical; it changes floating-point
summation order. Tests compare to direct pairs and a 90-digit Decimal oracle.
The implementation stores all distinct subset keys and bounds their count.

Interpretation limits
----------------------
The signed covariance estimator need not be positive semidefinite on a single
sample. Preserve and flag it; projecting it changes the estimator. Its
diagonal is mathematically nonnegative for these count variables, so a
negative diagonal indicates numerical cancellation. Zero observed motifs
can yield both zero estimate and zero estimated variance despite a positive
population count. Unbiased covariance does not certify Gaussian intervals.
The exact star counterexample in the validation directory demonstrates this.

Comparison with existing software
---------------------------------
R sampling::varHT evaluates the same variance when given unit values and
joint inclusion probabilities. R survey::ppscov accepts a sparse weighted
sampling-covariance matrix and svytotal returns joint covariance. These
capabilities are confirmed by actual execution, not merely documentation.
MotifSurvey's distinction is constructing the graph-specific calculation
without materializing all overlapping motif pairs, together with scoped input
checks, diagnostics and reproducible output. Sparse-pair methods can be faster
when the overlap matrix is sparse. The local vectorized-subset comparator is
another implementation of the same identity, not a published competing tool.

Linear contrasts and a full-population validation oracle
-------------------------------------------------------
For fixed coefficients a, the estimated total a^T T_hat is unbiased and has
unbiased variance estimate a^T V_hat a. For the closure example the vector
is (0, -0.25, 3). The resulting estimand is positive exactly when global
transitivity exceeds 0.25, provided that the population wedge count is positive.
This sign equivalence is not a significance test. A linear-contrast variance
estimate can be negative when the raw estimated covariance is indefinite,
even though the three count-variance diagonal entries are nonnegative.

For validation only, suppose the entire population graph is known. Let h(A,t)
count population type-t motif units containing nonempty vertex subset A.
The true design covariance is

  V = sum_A [product_(v in A) (1/p(v) - 1)] h(A) h(A)^T.

To derive this, the ordered-pair population summand is 1/pi(M intersection N)
minus one. Expand product_(v in M intersection N) [1 + (1/p(v)-1)] minus
one, then aggregate incidence counts by A. All coefficients are nonnegative.
This is a validation identity, not a new estimator. Computing it requires the
complete population, information that the sampled estimator does not use.

The example evaluates h for vertices and pairs using degree/common-neighbor
identities. An open-wedge triple has h=(0,1,0); a triangle triple has h=(0,3,1).
Tests compare the oracle with literal ordered pairs for all 64 four-node
graphs under two probability regimes and with an exhaustive sample distribution.
A separate weighted-degree/triple implementation also reproduces the recorded
sample contrast totals without using the subset covariance calculation.
