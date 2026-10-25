Release validation
==================

``compatibility_matrix.json`` records five executed Windows environments:
Python/NumPy 3.10.22/1.24.0, 3.12.14/1.26.4, 3.12.14/2.5.4,
3.14.9/2.5.4 and 3.15.1/2.5.4. Python 3.10 uses NetworkX 3.4.2;
the other cases use NetworkX 3.6.1. The Python 3.10, 3.14 and 3.15
runtimes were local Astral python-build-standalone installations managed by
uv 0.13.5. Their versions are recorded; no system Python was replaced.

For each interpreter, create a clean environment and run from a source copy::

    python -m venv .venv
    # Use that environment's Python for all subsequent commands.
    python -m pip install "numpy==VERSION" "path/to/wheel.whl[validation]"
    python -m pip check
    python -m unittest discover -s tests -v
    python -m motifsurvey examples/tiny.json
    python validation/followup/stress.py
    python validation/followup/validate_v020.py

Replace VERSION and the wheel path with the selected tested combination and
supplied wheel. Activation is shell-specific; invoking the environment's
Python by its full path is equivalent. The validation extra automatically
chooses the compatible NetworkX version. Install the source distribution in
another clean environment and repeat the unit suite and CLI example.

Each matrix case repeated the same 23 unit tests, 25 high-precision fixtures,
five expected rejections, 896 exact samples, 58 cached R-reference comparisons
and 60 public replays. Cached comparison means the R outputs were compared;
R was not rerun inside every Python environment. These repeated checks must
not be reported as additional independent scientific samples.

``public_study_v020.json`` contains all 600 email-network realizations replayed
under the current core. ``public_replay_summary.json`` compares all totals
and covariance entries to the preserved original 0.1.0 results at relative
and absolute tolerance 1e-8. Reproduce in a working source copy with::

    python validation/followup/public_study.py

The script verifies the downloaded data checksum before calculation and
writes ``validation/followup/public_study_results.json``. Preserve the supplied
historical record before rerunning. The 600 samples comprise 300 per design;
replaying their seeds does not increase the Monte Carlo sample size.

All performance claims in the paper use the original benchmark records and
their stated environment. The release replay ran concurrently with document
preparation and is used for numerical agreement, not new performance claims.
The CI action revisions are pinned in ``ci_action_versions.json``. No remote
workflow run, Linux test or macOS test is represented as complete.
