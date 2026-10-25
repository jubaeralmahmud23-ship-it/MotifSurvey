Planning precision for a network-closure survey
===============================================

Question
--------
At the same expected number of selected people, does emphasizing small or
large departments improve precision for network closure? The fixed reference
population is SNAP email-Eu-core projected to a simple undirected graph.
The complete 42-department roster supplies allocation information. This is
a design-planning experiment on an archived population; no claim is made
that the original email collection was a probability sample.

The estimand is D = 3 * triangles - 0.25 * wedges. For a population with at
least one wedge, its sign indicates whether global transitivity exceeds the
illustrative threshold 0.25. The estimator is linear and unbiased. Its variance
estimate is 9 V_triangle,triangle + 0.0625 V_wedge,wedge - 1.5 V_wedge,triangle.
Dropping the last term is not a valid calculation of contrast variance.
Individual raw estimates can be negative; they are retained, not clipped.
No normal approximation, confidence interval or hypothesis test is supplied.

Protocol and data
-----------------
``protocol.json`` was frozen internally before examining comparison results;
it is not an externally registered protocol. Three allocation scores are 1,
inverse square root of department size, and square root of department size.
Bounded calibration gives 201 expected selected vertices under each plan.
Vertex draws remain independent; this is not fixed-size stratified sampling.
All selected vertices, including isolates, and all induced edges are retained.

The pilot has 200 rounds; evaluation has 600 fresh rounds, with disjoint seeds.
Each round shares uniforms across the three plans, creating paired comparisons.
There are 2,400 plan realizations over 800 independent uniform vectors, not
2,400 mutually independent draws. The two pilot choices are frozen before
evaluating either choice. All three plans and their exact reference variances
are reported. Expected induced-edge costs are not equalized.

Source page: https://snap.stanford.edu/data/email-Eu-core.html
Both raw gzip downloads are verified against SHA256 before calculation.
Department labels SHA256:
e5abe5b4581a480032a63adcf2576c161785f45692642c6ebb0b1276f0f33669
Raw network and roster files are not redistributed or relicensed.

Reproduce
---------
From the extracted source directory::

    python -m pip install ".[validation]"
    python -m unittest discover -s tests -v
    python validation/design_case/run.py --cache data-cache --output fresh-design-results

Use a new output directory; the command refuses to overwrite previous results.
The checksum of the protocol and the executed core are recorded. Results include
the node probabilities, every estimate/covariance, pilot choices and evaluation
summary. Runtime fields are provenance only, not controlled performance claims.
``results/`` contains the executed release evidence. Bootstrap intervals measure
Monte Carlo uncertainty of paired MSE ratios, not uncertainty about D.

Inspect a previously computed MotifSurvey JSON result::

    python examples/closure_design.py result.json --threshold 0.25

A fixed follow-up, specified after seeing finite Monte Carlo discrepancies
in the original evaluation, checks true variance with a separate direct-total
formula. It reproduces all 2,400 existing contrasts before running 10,000 fresh
paired rounds (30,000 total-only calculations)::

    python validation/design_case/confirm_variance.py --cache data-cache --results fresh-design-results

These calculations use weighted-degree identities and explicit triangle triples,
not MotifSurvey's subset covariance algorithm. They are not additional covariance
estimator runs. The original evaluation and its discrepancies remain in place.
The archived follow-up MSE/exact-variance ratios are 1.073, 1.042 and 1.037;
their Monte Carlo errors are reported, without a perfect-calibration claim.

The exact-population oracle in the example has access to the entire reference
graph. Actual sampled estimation does not use its hidden population information.
The oracle uses specialized vertex/pair/triple incidence counts, independently
tested against literal motif pairs for 128 graph/probability combinations.
An exhaustive distribution over all 16 samples of a further graph checks
the contrast expectation, its true variance, and mean estimated variance.

Interpretation
--------------
The case demonstrates a reproducible precision comparison for one reference
network and three specified allocations. It does not establish an optimal
allocation, generalize to other networks, demonstrate external adoption, or
replace independent statistical review. The exact diagonal-only criterion
and its finite pilot can rank plans differently through Monte Carlo noise;
the results must not be represented as universal ranking reversal.
