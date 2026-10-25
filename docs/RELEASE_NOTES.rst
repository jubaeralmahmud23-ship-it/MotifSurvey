MotifSurvey 0.2.2
=================

This patch fixes probability input validation. NumPy Boolean and complex
values, byte strings and non-scalar arrays are rejected rather than converted
to probabilities. Overflowing integers now produce the documented ValueError
and a clean CLI error. Valid Python/NumPy real scalars, Decimal and Fraction
values remain supported. The closure-profile example applies corresponding
checks to thresholds.

Three regression tests raise the package suite to 31 tests; the four research
companion tests remain separate. Counting, weighting and covariance arithmetic
are unchanged for valid inputs. Historical studies and benchmarks retain their
original labels and are not relabeled as fresh experiments or timings.

This release includes the research companion, all 27 sensitivity settings,
native R interoperability records and the graph-to-survey workflow guide.
The frozen v0.2.1 release is preserved. The new release workflow builds and
checks its own wheel after the complete platform matrix passes.

Known limitations remain: independent vertex selection with known positive
probabilities, complete induced edges, simple undirected graphs, three motif
types, explicit enumeration budgets and no automatic confidence intervals.
No external adoption or independent human statistical review is claimed.
