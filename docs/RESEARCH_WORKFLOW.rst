Using MotifSurvey for a network research question
================================================

Who should use it
-----------------
Use MotifSurvey when the scientific target is an edge, wedge or triangle
total in a fixed population network, or a linear combination of those totals,
and the observed graph comes from documented independent vertex selection.
You must know each selected vertex's inclusion probability and observe every
edge between selected vertices. Include selected isolates. Sampling-design
assumptions come from the study design, not from the software.

One concrete question is whether a survey allocation gives useful precision
for network closure at an affordable expected number of selected vertices.
An archived network can serve as a reference population for planning. An
actual sampled analysis needs only the sample and its probabilities; the
exact population oracle used to evaluate this example is not available from
a sample alone.

What the package replaces
-------------------------
For this same question, a graph-to-R-survey workflow requires the researcher to:

1. Enumerate edges, centered wedges and triangles with consistent multiplicity.
   A triangle supplies three wedges, even though their vertex sets coincide.
2. Turn those motifs into survey units and compute each inclusion probability
   as the product of its vertices' probabilities.
3. Compute overlap-dependent joint inclusion information. Independent vertex
   draws do not imply independent motif observations.
4. Supply the resulting weighted overlap covariance to survey::ppscov and
   survey::svydesign, then obtain totals and joint covariance with svytotal.
5. Apply the desired contrast and retain the inputs, software version and
   diagnostics needed to audit the calculation.

Our executed native adapter is validation/followup/native_compare.R. It makes
the alternative reproducible and shows that R survey can produce the same
answer. MotifSurvey supplies these motif-specific steps behind a graph API,
using subset aggregation instead of constructing a matrix of overlapping
motif pairs. The implementation adds design declarations, input/resource
checks, canonical input hashes and explicit warnings about unstable inference.
These are concrete software features; no measured human-time savings or
unique statistical capability is claimed.

Choosing between tools
----------------------
* Use NetworkX or igraph for counts and structure of the graph you already
  observe, particularly when no population sampling inference is intended.
* Use MotifSurvey for this supported sampling design when you need the three
  weighted totals, their joint covariance and a checked graph input/output
  workflow. It still enumerates motifs and can exceed its resource budgets.
* Use R survey with a suitable graph adapter when working in an established
  R survey pipeline or needing its broader survey-analysis facilities. Other
  sampling designs require their own correct motif inclusion information;
  MotifSurvey does not supply that extension.
* A local vectorized implementation of the same subset identity is often
  faster. Choose it for kernel experimentation if you will supply and maintain
  the surrounding checks yourself. It is a comparison implementation, not an
  independently released competitor.

From a sample to a reusable result
---------------------------------
Install the released 0.2.1 wheel and obtain this research companion source.
With an input conforming to docs/API.rst, run from the companion directory::

    python -m motifsurvey input.json --output result.json
    python examples/closure_profile.py result.json --thresholds 0.15 0.25 0.35 --output profile.json

The first command enumerates motifs once. The second reuses its three totals
and full covariance for D(t) = 3 T_triangle - t T_wedge. For a population with
a positive wedge total, the sign of D(t) indicates whether transitivity
exceeds t. The estimated variance includes the covariance term
-6 t V_wedge,triangle. No second graph enumeration or new observations are
needed for these linear combinations.

Specify scientifically meaningful thresholds before interpreting results.
The supplied grid is an exploratory sensitivity example; choosing thresholds
after seeing estimates does not produce valid selection-adjusted inference.
The profile retains negative variance estimates, the input hash and source
warnings. It supplies no confidence intervals, significance tests, coverage
or power guarantees. Do not interpret a positive point estimate as evidence
that the population threshold has been established.

What the executed example teaches
---------------------------------
At 201 expected selected vertices, replacing individual probabilities by their
mean gives exact triangle-total biases of +26.98% for small-department emphasis
and -7.12% for large-department emphasis. This is a diagnostic misuse, not a
competing package or a recommended estimator. Correct HT weighting is unbiased
under the assumed design, although it need not minimize mean squared error:
the biased diagnostic has lower empirical closure-contrast MSE in one of the
three archived allocations. Every result is retained in the companion.

Exact standard deviations of D(0.25) under uniform selection are 68,186,
37,846 and 20,851 at expected vertex budgets 100.5, 201 and 402. The reference
target is only 20,579. Thus this planning exercise exposes substantial sampling
uncertainty even at the largest budget; it does not establish an adequately
powered design. Uniform selection is the most precise of the three candidates
for this target, but at threshold 0.35 and budget 402, small-department emphasis
has 15.0% lower exact variance. Allocation preference depends on the estimand.
These are exact calculations on one reference network, not general optimality.

The first evaluation sample under each allocation was also analyzed by the
unmodified native R survey function through the supplied graph adapter. Three
totals, their complete covariance, and all nine contrast/variance pairs agreed
within the specified 1e-8 tolerance (maximum scaled profile error 2.76e-13).
This comparison establishes interoperability, not statistical superiority.

See validation/research_usefulness/README.rst for the complete protocol, all
27 sensitivity settings, archived-sample reanalysis and reproduction commands.
No external adoption, user study or independent human statistical review is
claimed. The companion adds examples and evidence around the frozen 0.2.1
software; it does not change its installed estimator or replace its release.
