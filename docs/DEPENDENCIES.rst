Dependencies and data
=====================

Runtime
-------
The package imports Python's standard library and NumPy (declared >=1.24).
NumPy is distributed separately; its BSD license and notices for bundled
components apply to its wheels. MotifSurvey contains no vendored NumPy code
or binary. The package itself is distributed under the MIT License in
``LICENSE`` and ``Licence.txt``; these files are identical.

Validation only
---------------
NetworkX supplies independent motif counts and graph fixtures. The validation
extra selects NetworkX 3.4.2 on Python 3.10 and 3.6.1 on Python >=3.11, to
respect their declared Python requirements. NetworkX uses the BSD license.
The performance adapters also use SciPy 1.18.1 (BSD license); Windows native
process measurements use psutil 7.2.2 (BSD license). These are separate
installations, not runtime requirements of MotifSurvey.

Native comparison scripts invoke separately installed R 4.6.1, sampling 2.11
(GPL >=2), survey 4.5 (GPL-2 or GPL-3), Matrix 1.7-5 and jsonlite 2.0.0.
Their upstream licenses apply to those installations. The source archive
contains our adapter scripts and outputs, not the R packages or their source.
Comparing against a package is not a claim to relicense it.

Official sources, checked 25 October 2026:

* https://numpy.org/doc/stable/license.html
* https://networkx.org/documentation/stable/reference/legal.html
* https://scipy.org/about/
* https://github.com/giampaolo/psutil/blob/master/LICENSE
* https://cran.r-project.org/package=sampling
* https://cran.r-project.org/package=survey

Public datasets
---------------
The experiments download from SNAP; no raw network file is redistributed.
Data availability does not imply that MotifSurvey's MIT License covers the
data. Consult the source pages and acknowledge the original researchers.
Stored node IDs are used only to specify reproducible synthetic probabilities.

email-Eu-core:
https://snap.stanford.edu/data/email-Eu-core.html
SHA256 of the downloaded gzip:
4b47acdb80197b085fe63c819c357ae488131ee904ed93d1b219a68b0f9e245f
Remove self-loops, merge reciprocal records and retain all endpoints before
counting. This produces a simple undirected artificial population with 1,005
vertices and 16,064 edges. Cite Yin, Benson, Leskovec and Gleich (2017), Local
Higher-order Graph Clustering, and the source page's additional citation to
Leskovec, Kleinberg and Faloutsos (2007), Graph Evolution: Densification and
Shrinking Diameters.

The closure-design case also downloads the department roster from the same
email-Eu-core page. Its gzip SHA256 is:
e5abe5b4581a480032a63adcf2576c161785f45692642c6ebb0b1276f0f33669
The roster is used for allocation scores; raw labels are not bundled.

Earlier ca-GrQc experiment:
https://snap.stanford.edu/data/ca-GrQc.html
SHA256 of the downloaded gzip:
a254442cdf5d684712578b630c2e0d7543518ab154ef2341cabb607572ce7230
Remove loops and merge repeated undirected edges, retaining endpoints.
This produces 5,242 vertices and 14,484 edges. Cite the source page and
Leskovec, Kleinberg and Faloutsos (2007).

No inference about either original data-collection mechanism is made.
