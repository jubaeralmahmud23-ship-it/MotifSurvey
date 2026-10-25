Release checks
==============

Version 0.2.2 fixes invalid scalar conversions. Its local installed wheel
passed 31 core tests, four companion tests and the complete release checker.
The exact-fraction audit covers 128 sampling laws and 2,048 sample calculations;
all 2,400 archived contrast records were recomputed. See validation/audit/.
All 13 platform checks and the publication job passed for tag v0.2.2:
https://github.com/jubaeralmahmud23-ship-it/MotifSurvey/actions/runs/40192521004
See HOSTED_CHECKS_v022.json and PUBLIC_RELEASE_v022.json for downloaded
reports, immutable source identity and verified public-asset checksums.
The new core SHA256 is:
b26c1e4531c6b96e9bfc3b4c34ec55bd78e8cd032498b2cec25fb2f0c31c4cc1

The author confirmed review, contribution roles, sole authorship and no
competing interests on 25 October 2026. See MANUSCRIPT_STATUS.rst for the
current manuscript status, including confirmed unpublished status and no
simultaneous submission.

Historical records below retain their original versions.

The release preparation report in ``validation/release/`` records executed
installation and numerical checks. The historical 0.2.1 core hash is:
4412fcb83e74aaf492bd3eaf1ea7509a5ec51d4d31e957ebb6dc3cf177d51c31
Compared byte-for-byte after normalizing line endings with the archived 0.2.0
wheel, the only source changes import and report the shared version string.
No estimation arithmetic changed. The old core hash is retained in historical
0.2.0 records, not asserted to describe the current file.

The earlier follow-up records cover 23 unit tests, 896 exhaustive samples,
58 distinct native graph/tool comparisons (78 runs), 25 high-precision cases,
five deliberate resource/numeric rejections, 60 public replays and the large
star check. Their timing environments and original version labels are retained.

All 13 hosted jobs passed on 25 October 2026: Windows, Linux and macOS with
Python 3.10, 3.12, 3.14 and 3.15, plus the minimum NumPy 1.24.0 job.
Each report records 28 unit tests, high-precision/exact checks and cached
replays. See HOSTED_CHECKS_v021.json for actual environments and job links.
These compatibility repetitions are not new independent scientific experiments.

Before submitting the manuscript
--------------------------------
1. Confirm the author names, order, affiliations, contact, funding, competing
   interests, contribution roles and copyright ownership.
2. Review the code, derivation, evidence and manuscript. Independent statistical
   review is recommended, and is not claimed to have occurred.
3. Public source, the versioned release and the 13-job CI matrix are complete.
   An archival deposit may be added where appropriate. No DOI is claimed.
4. Recheck the current journal guide and submission portal. The official
   template retrieved on 25 October 2026 specifies five main sections,
   GitHub, README.md, Licence.txt, 4,000 words and at most six figures.
   The full ScienceDirect guide returned HTTP 403. Fees and any additional
   portal requirements have not been verified.
5. Confirm all declarations and final manuscript approval; do not submit a
   draft containing pending metadata or unconfirmed author statements.

These steps do not imply that acceptance, external adoption or all-platform
support has already been established. The software is a bounded research
implementation for a known sampling design.

Historical 0.2.0 compatibility matrix
-------------------------------
Python/NumPy: 3.10.22/1.24.0, 3.12.14/1.26.4, 3.12.14/2.5.4,
3.14.9/2.5.4 and 3.15.1/2.5.4. All passed the same 23 unit tests,
25 Decimal cases, five expected rejections, 896 exact samples, 58 cached
native comparisons and 60 selected public replays. A separate full replay
compared all 600 current-version samples to the archived original results.
See validation/release/compatibility_matrix.json and public_replay_summary.json.

Executed 0.2.1 compatibility matrix
----------------------------------
The same five Python/NumPy combinations passed the installed-wheel checker
with all 28 unit tests, 25 high-precision cases, five expected rejections,
896 exact samples, 58 cached native comparisons, 60 public replays and pip
dependency checks. Source installation also passed the 28-test suite. See
``validation/release/compatibility_matrix_v021.json``. Validation fixes BLAS
thread counts to one to avoid oversized thread pools on small machines.
An earlier unrestricted local subprocess failed to allocate a BLAS thread
pool during concurrent work; no numerical test failure was suppressed.

Workflow syntax checking with actionlint 1.7.12 passed. Actual hosted job
results are available at:
https://github.com/jubaeralmahmud23-ship-it/MotifSurvey/actions/runs/40166324430

The initial combined branch/tag push triggered the package matrix but did not
start the separate Versioned release workflow. The initial release was
published through the API after verifying every job and report. The attached
local artifacts match their verified Windows release checks; hosted runners
independently built the same tagged source. No automated publication-job
execution or byte-identical hosted installation of the attached wheel is claimed.

Release: https://github.com/jubaeralmahmud23-ship-it/MotifSurvey/releases/tag/v0.2.1
Immutable commit: e5235add0fb3cb8fe768e6baf846c894438c7b56
The tag remains unchanged; this branch records post-deposition evidence.
