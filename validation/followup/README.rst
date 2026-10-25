Extended validation and reproduction
====================================

This directory preserves the extended-validation decision and records.
Version 0.2.0 is the delivered core. The 600 public thinnings originally used
installed 0.1.0; the follow-up replayed 60 samples with 0.2.0 and passed a
fresh 896-sample exact oracle. Subsequent release preparation replayed all
600 samples and tested five Python/NumPy environments: see ../release/.
Repeated replays are not additional independent Monte Carlo samples.

Install from the extracted package root::

    python -m pip install ".[validation]"
    python -m pip install scipy==1.18.1 psutil==7.2.2
    python -m unittest discover -s tests -v
    python validation/followup/stress.py
    python validation/followup/public_study.py
    python validation/followup/validate_v020.py

The public study downloads SNAP email-Eu-core, hashes it, removes loops and
projects it to an undirected simple graph. It writes local JSON results and
performs 600 synthetic thinnings. It validates the stated artificial sampling
design, not the source's data collection. The replay script can also fetch the
dataset directly. In a fresh rerun, public-study results identify the installed
package version. Archived original outputs are evidence from the stated run,
not promises that timings or last floating-point bits will repeat.

Windows performance experiments::

    python validation/followup/benchmark.py
    python validation/followup/paired_versions.py
    python validation/followup/large_star.py

``benchmark.py`` measures process peak working set using Windows APIs and is
not a cross-platform memory harness. The statistical code itself does not
use these APIs. ``paired_versions.py`` compares the delivered implementation
with our frozen original 0.1.0 core, in randomized order with seven repeats.
The frozen core is not a statistical oracle. ``large_star.py`` deliberately
raises the resource budget and checks an independent closed-form expression.
No huge quadratic comparator is run for that star.

Native R comparisons
---------------------
Install R and make Rscript available on PATH. The actual benchmark used R
4.6.1 with sampling 2.11, survey 4.5, Matrix 1.7-5 and jsonlite 2.0.0. The
installer is not distributed here. ``runtime_fetch_log.json`` records its
official URL and checksums. Future package installation may resolve newer
versions; record and pin them for an exact-version replication.

From the package root::

    Rscript --vanilla validation/followup/install_r_packages.R validation/followup/r-library
    python validation/followup/prepare_native.py
    Rscript --vanilla validation/followup/native_compare.R validation/followup/r-library validation/followup/native_fixtures.json validation/followup/native_results.json
    python validation/followup/check_native.py

These execute unmodified installed functions through our graph-to-motif
adapter. Timing distinguishes enumeration, setup and native inference. R's
clock resolution may give zero for very short stages. The native ``survey``
design API rejected the first one-unit fixture; the comparison fixtures were
changed to have at least two edges, and all experiments rerun. The package's
single-edge behavior is independently unit tested.

On Windows, ``native_memory.py`` measures the native process tree. Set the
MOTIFSURVEY_RSCRIPT environment variable to the Rscript executable if it is
not on PATH; set MOTIFSURVEY_RLIB to the installed package library when using
a different location. It reads peak working set for the launcher and its
descendants, reporting the largest process peak. The initial launcher-only
measurement was invalid and is excluded. Full native peaks include libraries;
they are not an isolated measure of the covariance kernel's memory.

Scope of the supplied evidence
------------------------------
* Native R comparisons: 78 runs, 30 graphs, 58 distinct graph/tool pairs.
* Seven stronger Python workloads. Two potentially huge sparse-pair matrices
  were deliberately not run and have no measured speedup reported.
* 25 high-precision cases and five expected rejection cases.
* 600 second-public-graph thinnings, plus the earlier 60 ca-GrQc thinnings.
* 23 current unit tests, 896 current exact samples and 60 current public replays.
* A closed-form 499,500-unit star and paired version timings.

No human statistical review, user study, public release, another-platform run
or confidence-interval coverage validation is claimed. ``THEORY.rst`` is an
internal derivation audit. Official journal-guide access was still blocked.
The software-selection score is not an acceptance probability.
