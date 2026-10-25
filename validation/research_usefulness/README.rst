Research usefulness companion to MotifSurvey 0.2.1
==================================================

This supplement adds an executable research workflow and a sensitivity
analysis to the previously released package. The installed src/motifsurvey
files and the v0.2.1 tag are unchanged. It is a retrospective expansion of
the same reference-network case, not independent empirical replication.

Read docs/RESEARCH_WORKFLOW.rst for the research question, tool-selection
guidance and graph-to-survey conversion. protocol.json specifies the full
grid before this supplementary calculation: three expected vertex budgets,
three allocation plans and three closure thresholds. All 27 settings appear
in results/summary.json. The complete graph supplies exact variance and
expectation oracles; no hidden population information enters the sample
estimator. Bounds of 0.05 and 0.8 can change the allocations as budget changes.
Costs are equalized in expected vertices, not edges or realized sample size.

The diagnostic mean-probability correction divides observed counts by the
second or third power of the mean vertex probability. It is correct for
uniform sampling; it is generally biased for unequal probabilities. It is
not an existing software baseline. Exact bias and all 600 archived evaluation
realizations per plan are reported. The reanalysis uses no new random draws.
Unbiasedness does not guarantee lower MSE; the diagnostic has lower empirical
closure-contrast MSE under large-department emphasis, despite its bias.

The native comparison uses the first evaluation seed, 2026103500, under all
three plans at budget 201, not samples selected for favorable results.
Unmodified R survey 4.5, Matrix 1.7-5, sampling 2.11 and jsonlite 2.0.0 were
loaded under R 4.6.1. Only survey is used in these three runs. The original
native_compare.R adapter is unchanged; its hash is recorded. Its measured
single-run stage times are provenance, not new performance comparisons.
The nine derived contrast/variance pairs are postprocessed from each tool's
three totals and covariance. Cached comparisons in CI are replays, not R runs.

Reproduce from this companion's root directory
---------------------------------------------
Install the versioned 0.2.1 wheel and its validation dependencies, then::

    python -m pip install "networkx==3.6.1"
    python validation/research_usefulness/run.py --cache data-cache --output fresh-usefulness
    Rscript validation/followup/native_compare.R /path/to/R/library fresh-usefulness/native_fixtures.json fresh-usefulness/native_results.json
    python validation/research_usefulness/check_native.py fresh-usefulness
    python -m unittest discover -s validation/research_usefulness -p "test_*.py" -v

Python 3.10 uses NetworkX 3.4.2 instead of 3.6.1. R's library directory must
contain the four packages above. The runner refuses an existing output
directory. Raw datasets are fetched from SNAP and checksum verified using
the original closure_design.py loader. Native graph fixtures are regenerated
from those data and the fixed seed; raw graph/roster downloads are not bundled.

Four companion tests check exact first moments against all 16 samples of all
64 four-vertex graphs under two unequal-probability regimes (2,048 sample
calculations); retain negative profile variances; reject invalid order and
thresholds; and replay the archived native comparisons. They are additional
to the original package's 28 tests. execution_notes.json records development
repairs to this new diagnostic; the released estimator was unaffected.

The profile is an example script, not an added installed API. It reuses the
JSON result and does not estimate a ratio, construct intervals or adjust for
post-selection. Profile thresholds and budget comparisons are sensitivity
analyses on the same archived reference graph. No external users, validation
of original data collection, or independent human review is asserted.
