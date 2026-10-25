Reference and claim audit, 25 October 2026
=========================================

Nine manuscript bibliography entries were checked against primary sources.
This audit distinguishes bibliographic verification, full-text inspection,
documentation inspection and actual execution. It does not claim every cited
paper was independently replicated or that all possible competitors were found.

1. Horvitz and Thompson (1952), JASA 47(260), 663-685,
   DOI 10.1080/01621459.1952.10483446. Publisher metadata and abstract verified;
   publisher full text was inaccessible. The formula is independently derived
   in THEORY.rst and checked against the R sampling manual and executable.
   https://www.tandfonline.com/doi/abs/10.1080/01621459.1952.10483446

2. Klusowski and Wu (2018), Counting Motifs with Graph Sampling, PMLR 75,
   1966-2011. Official proceedings metadata and abstract verified. It studies
   induced motif counts; it supports network-sampling background, not a claim
   that its estimand is identical to our non-induced centered wedges.
   https://proceedings.mlr.press/v75/klusowski18a.html

3. Bhattacharya, Das and Mukherjee, Annals of Statistics 50 (2022), 987-1011,
   DOI 10.1214/21-AOS2134. The consulted formula is explicitly attributed to
   arXiv v2 dated 21 December 2024. Proposition 2.4 estimates the variance of
   the unscaled sampled motif count. Dividing by p^(2k), for k motif vertices,
   gives the equal-probability HT overlap expression used here. Normality in
   that paper requires additional conditions; we do not assert they hold in
   every package input or provide automatic confidence intervals.
   https://arxiv.org/html/2011.03026v2

4. NetworkX: executed version 3.6.1; citation now pins the 3.6 documentation
   series instead of a moving stable URL. triangles counts observed-graph
   triangles and does not itself correct independent population subsampling.
   https://networkx.org/documentation/networkx-3.6/reference/algorithms/generated/networkx.algorithms.cluster.triangles.html

5. python-igraph: documentation inspected at version 0.11.8, not executed.
   motifs_randesu_estimate extrapolates motif totals from sampled starting
   vertices and randomized search within the supplied graph. This is different
   from inference about an unobserved population from a pre-observed induced
   vertex sample. No claim that igraph never samples vertices is made.
   https://python.igraph.org/en/0.11.8/api/igraph.GraphBase.html#motifs_randesu_estimate

6. Tille and Matei, sampling 2.11 (2025). The official manual's varHT method=1
   is the general unbiased HT variance. Method=2 is the fixed-size
   Sen-Yates-Grundy form and is not the comparison used. Actual native runs
   and distinct graph/tool comparison counts are archived separately.
   https://cran.r-project.org/web/packages/sampling/sampling.pdf

7. R survey 4.5 (2026). Official CRAN metadata and ppscov documentation checked.
   weighted=TRUE supplies sampling covariance divided by joint inclusion
   probabilities. svytotal reproduces the joint result after our graph-to-motif
   conversion. Finite agreement supports numerical correctness, not a universal
   proof. Actual native runs are distinct from cached CI replays.
   https://cran.r-project.org/web/packages/survey/index.html
   https://github.com/cran/survey/blob/master/man/HR.Rd

8. Mandal and Chatterjee, arXiv:2407.19191v2 (13 April 2026), title and authors
   verified. This is model-and-design inference in a stochastic block model.
   The comparison is limited to selected author source previously inspected at
   d72ae36050fa5715523a7abc31cbb2c14f3e5187, not execution of the full repository.
   https://arxiv.org/abs/2407.19191v2
   https://github.com/AnirbanM13/subgraph-count-estimate/tree/d72ae36050fa5715523a7abc31cbb2c14f3e5187

9. SNAP email-Eu-core: official dataset page and attribution checked. The input
   is directed; removing loops and merging reciprocal edges defines a different,
   undirected reference population. Controlled thinning tests our declared
   sampling design and does not validate the original email collection design.
   The source and department-label hashes are retained in the study records.
   https://snap.stanford.edu/data/email-Eu-core.html

No bibliography entry was found to be fabricated. Documentation versions and
the variance rescaling attribution were clarified. Exhaustive reference
verification beyond accessible primary sources is not claimed.
