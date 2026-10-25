# MotifSurvey

Estimate population edge, wedge and triangle totals from an induced graph obtained by independently sampling vertices with known inclusion probabilities. MotifSurvey returns the three Horvitz-Thompson totals and their joint estimated design covariance. Probabilities may differ between vertices.

Version **0.2.2** fixes invalid probability conversions found during the mathematics and code audit. Valid-input estimation arithmetic is unchanged. The frozen [0.2.1 release](https://github.com/jubaeralmahmud23-ship-it/MotifSurvey/releases/tag/v0.2.1) and its evidence remain available. The [0.2.2 release workflow](https://github.com/jubaeralmahmud23-ship-it/MotifSurvey/actions/workflows/release.yml) runs the full platform matrix before publishing assets. No archival DOI has been assigned.

## Install

Use Python 3.10 or later:

```sh
python -m pip install .
python -m motifsurvey examples/tiny.json
```

From a supplied wheel, use `python -m pip install path/to/motifsurvey_prototype-0.2.2-py3-none-any.whl`. The distribution name remains `motifsurvey-prototype`; the import and command name are `motifsurvey`. NumPy is the only runtime dependency. No compiled extension is included.

## A small example

```python
from motifsurvey import estimate

result = estimate(
    nodes=["a", "b"],
    edges=[("a", "b")],
    probabilities={"a": 0.5, "b": 0.5},
    design="independent_node",
    complete_induced=True,
)
print(result["estimated_totals"])      # [4.0, 0.0, 0.0]
print(result["estimated_covariance"][0][0])  # 12.0
```

The edge has observation probability 0.25, so its contribution to the population estimate is four. This is an estimate from a sample, not a statement that the population contains exactly four edges.

For JSON input, the fields are the same as the function arguments. `python -m motifsurvey input.json --output result.json` writes a new file and refuses to overwrite an existing one. The CLI exits with status 2 for rejected input. See [API and input format](docs/API.rst) for output fields and resource budgets.

## What must be true about the sample

The population is a fixed finite simple undirected graph. Every population vertex has a known, positive selection probability; vertices are selected independently; and all edges between selected vertices are observed. The input lists all selected vertices, including isolates, and exactly one probability per selected vertex. Vertex IDs are nonempty strings. Duplicate edges and self-loops are rejected.

These are sampling-design assumptions supplied by the researcher. They cannot be verified from the observed graph. Convenience samples, estimated probabilities, fixed-size samples, ego networks, snowball samples and graphs with missing induced edges are outside the supported scope.

Counts are non-induced. A wedge is a center with an unordered pair of neighbors; a triangle contains three wedges. Zero observed motifs do not establish population absence. A raw estimated covariance can be indefinite, and the output flags this without changing the estimator. The software does not provide automatic confidence intervals or certify a normal approximation.

## Computation and comparisons

For a research-oriented walkthrough, start with [Choosing tools and analyzing a sampled network](docs/RESEARCH_WORKFLOW.rst). It explains the graph-to-survey conversion, when R survey or a graph counter is the better fit, and how to reuse one result for closure contrasts. The [research companion](validation/research_usefulness/README.rst) adds all 27 budget/allocation/threshold settings and actual R comparisons on the same sampled graphs. These additions leave the released 0.2.1 estimator unchanged.

The [pinned companion](https://github.com/jubaeralmahmud23-ship-it/MotifSurvey/tree/e15eecff5f70a1616c53816bb75fdf7c2acffb4a) passed a [separate 13-job matrix](https://github.com/jubaeralmahmud23-ship-it/MotifSurvey/actions/runs/40170607039), including its four additional tests. [Downloaded environment reports](docs/RESEARCH_COMPANION_CHECKS.json) preserve the evidence. These tests replay cached native outputs; the three actual R runs are recorded separately in the study.

The estimators are established survey-sampling theory. MotifSurvey aggregates weighted motif incidences by nonempty vertex subsets and computes the covariance in blocks of at most 512 keys. This avoids a motif-by-motif covariance matrix while still enumerating every motif. It is an independent implementation, without a claim of a new estimator or first use of inclusion-exclusion.

Default limits are 200,000 motif instances and 1,000,000 subset keys. These are count budgets, not a guaranteed RAM bound. Dense graphs and high-degree vertices can exceed them. A graph that is too large is rejected rather than silently sampled or truncated.

Native comparisons executed `sampling::varHT` (2.11) and `survey::svytotal` with `ppscov` (4.5) in R 4.6.1. These packages implement the same underlying design-based covariance after graph-to-motif conversion. Our adapters and their setup costs are supplied separately. In 78 comparisons across 30 fixtures, results agreed within a tolerance of 1e-8. The R packages are not included in this distribution.

Strong Python baselines include a local SciPy sparse-overlap adapter and a local adapter using the same subset identity. Sparse overlaps are faster on rings and one random graph; MotifSurvey saves time and peak process memory when many motifs overlap. The same-identity adapter is often faster but lacks the package's input checks and provenance output. The benchmark records report these differences rather than claiming a universal speedup.

## Reproduce the checks

The [closure-design case](validation/design_case/README.rst) asks which independent sampling plan estimates network closure most precisely at an equal expected vertex budget. It uses department-roster information, a separate pilot and evaluation, and all three candidates. The linear contrast uses the returned joint covariance; it is not a confidence interval or a general missing-network correction.

```sh
python -m pip install ".[validation]"
python -m unittest discover -s tests -v
python validation/followup/stress.py
python validation/followup/validate_v020.py
python tools/check_release.py --output check-results
```

Run these commands from the extracted source directory. The replay command downloads and verifies a small public SNAP dataset if it is absent. Scripts write result JSON next to themselves: use a working copy if preserving the supplied evidence. [Reproduction instructions](validation/followup/README.rst) describe native R comparisons, public thinning experiments and platform-specific performance measurements. [Theory](THEORY.rst) derives the estimator. [Release checks](docs/RELEASE_CHECKS.rst) records local and hosted checks and identifies the remaining human review.

The public experiments treat the stored networks as fixed artificial populations and apply controlled vertex sampling. They do not validate the sampling mechanism of the original data collection. All 600 email-network samples were replayed under 0.2.0 during release preparation and agreed with the archived 0.1.0 calculations. These are the same controlled samples, not 600 new independent realizations. Five fresh Windows environments also passed the unit, exact-sampling, high-precision and selected-replay checks.

## License, attribution and support

The package and supplied original validation scripts use the [MIT License](Licence.txt). Third-party libraries and public datasets retain their own terms; see [dependencies and data](docs/DEPENDENCIES.rst). Datasets and third-party package source are not bundled.

Maintainer: Jubaer Al Mahmud, Independent researcher, Dhaka, Bangladesh; jubaeralmahmud23@gmail.com. The [repository](https://github.com/jubaeralmahmud23-ship-it/MotifSurvey) and [versioned release](https://github.com/jubaeralmahmud23-ship-it/MotifSurvey/releases/tag/v0.2.1) are public and their downloads have been verified. Citation metadata is in `CITATION.cff`; no archival DOI has been assigned. No external user adoption or independent human statistical review is claimed.

Primary references:

- Horvitz and Thompson (1952), [A generalization of sampling without replacement from a finite universe](https://doi.org/10.1080/01621459.1952.10483446).
- Klusowski and Wu (2018), [Counting Motifs with Graph Sampling](https://proceedings.mlr.press/v75/klusowski18a.html).
- Bhattacharya, Das and Mukherjee, [Motif Estimation via Subgraph Sampling: The Fourth Moment Phenomenon](https://arxiv.org/abs/2011.03026v2), version 2, 2024. The uniform-probability overlap variance appears in Proposition 2.4, with HT rescaling by p^(-2k) for k-vertex motifs. This version is identified explicitly because the printed expression differs from an earlier published-copy PDF.
- [R sampling](https://cran.r-project.org/package=sampling) and [R survey](https://cran.r-project.org/package=survey).

## Version 0.2.2 audit

The [audit](validation/audit/README.rst) checks 2,048 sample calculations using exact fractions and recomputes all 2,400 archived contrast records. The local installed package passed 31 core tests, four companion tests, the high-precision and exhaustive checks, and cached replays. [Reference checks](docs/REFERENCE_AUDIT.rst) explain the sources and scope of each comparison. Historical studies retain their original versions and were not rerun as new experiments.

## Historical 0.2.1 release provenance

The tested source is commit `e5235add0fb3cb8fe768e6baf846c894438c7b56` (tag `v0.2.1`). [Hosted reports](docs/HOSTED_CHECKS_v021.json) record the actual environments and numerical checks. The [release record](docs/PUBLIC_RELEASE_v021.json) records anonymously verified downloads and hashes. Hosted jobs built the tagged source; the attached local wheel was separately verified on Windows. Compatibility replays are not new scientific samples.

The tag and downloadable distributions retain the documentation as it stood at source freeze. This branch adds the subsequent publication and hosted-test record without changing the tagged software. The initial combined branch/tag push started the package matrix only; publication followed explicit verification of its 13 successful jobs. The automated publication job was not executed for this initial release.
