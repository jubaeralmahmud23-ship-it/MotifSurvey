"""Classical HT estimation, evaluated through subset aggregation.

For observed motif units M,N, the covariance summand is
w_M w_N (1 - product(p_v for v in M intersection N)). Expand this
intersection term by inclusion-exclusion and aggregate weighted motif
incidences by their nonempty vertex subsets. This avoids motif-pair matrices.
The formula is established sampling theory; this is an independent implementation.
"""
from itertools import combinations
from math import prod, isfinite
import hashlib
import json
import numpy as np
from ._version import __version__

NAMES = ("edges", "wedges", "triangles")

class ResourceLimitError(ValueError):
    """The requested calculation exceeded an explicit resource budget."""

def estimate(nodes, edges, probabilities, *, design, complete_induced,
             max_instances=200_000, max_subsets=1_000_000):
    """Return HT totals and an unbiased estimated design covariance.

    nodes: unique string IDs for every observed vertex, including isolates.
    edges: unique unordered pairs; loops and parallel records are rejected.
    probabilities: positive known inclusion probabilities for exactly nodes.
    design must equal 'independent_node'; complete_induced must be True.
    These declarations cannot be verified from the observed graph. All
    population nodes must have positive inclusion probability by design.
    Wedges are non-induced, centered unordered neighbor pairs: a triangle
    contributes three wedges and one triangle. No normal intervals are returned.
    """
    if design != "independent_node" or complete_induced is not True:
        raise ValueError("Requires independent_node sampling and all induced edges observed")
    if any(isinstance(x,bool) or not isinstance(x,int) or x<1 for x in (max_instances,max_subsets)):
        raise ValueError("Resource limits must be positive integers")
    nodes=list(nodes)
    if any(not isinstance(v,str) or not v for v in nodes) or len(set(nodes))!=len(nodes):
        raise ValueError("Nodes must be unique nonempty string IDs")
    if set(probabilities)!=set(nodes):
        raise ValueError("Probability keys must match observed nodes exactly")
    p={}
    for v in sorted(nodes):
        value=probabilities[v]
        if isinstance(value,(bool,str)):
            raise ValueError("Probabilities must be numeric, finite and in (0,1]")
        try:value=float(value)
        except (TypeError,ValueError):raise ValueError("Invalid probability") from None
        if not isfinite(value) or not 0<value<=1:raise ValueError("Probabilities must be in (0,1]")
        p[v]=value
    adj={v:set() for v in sorted(nodes)};unique=set()
    for edge in edges:
        if isinstance(edge,(str,bytes)):raise ValueError("Each edge must be a pair")
        edge=tuple(edge)
        if len(edge)!=2:raise ValueError("Each edge must be a pair")
        u,v=edge
        if u not in adj or v not in adj or u==v:raise ValueError("Unknown endpoint or self-loop")
        key=tuple(sorted((u,v)))
        if key in unique:raise ValueError("Duplicate undirected edge")
        unique.add(key);adj[u].add(v);adj[v].add(u)
    # This bound is available before enumerating wedges or triangles.
    lower_bound=len(unique)+sum(len(a)*(len(a)-1)//2 for a in adj.values())
    if lower_bound>max_instances:raise ResourceLimitError("Edge/wedge count exceeds max_instances")
    totals=np.zeros(3);counts=np.zeros(3,dtype=np.int64);sums={}
    def add(vertices,t):
        vertices=tuple(sorted(vertices));pi=prod(p[v] for v in vertices)
        if pi==0 or not isfinite(1/pi):raise ValueError("Inclusion probability underflows or inverse overflows")
        weight=1/pi;counts[t]+=1
        if int(counts.sum())>max_instances:raise ResourceLimitError("Motif count exceeds max_instances")
        totals[t]+=weight
        for k in range(1,len(vertices)+1):
            for subset in combinations(vertices,k):
                if subset not in sums:
                    if len(sums)>=max_subsets:raise ResourceLimitError("Subset count exceeds max_subsets")
                    sums[subset]=np.zeros(3)
                sums[subset][t]+=weight
    for u,v in sorted(unique):add((u,v),0)
    for u in sorted(adj):
        for v,w in combinations(sorted(adj[u]),2):add((u,v,w),1)
    for u in sorted(adj):
        for v in sorted(adj[u]):
            if v>u:
                for w in sorted(adj[u]&adj[v]):
                    if w>v:add((u,v,w),2)
    covariance=np.zeros((3,3));absolute_terms=np.zeros((3,3))
    with np.errstate(over='raise',invalid='raise'):
        try:
            # Bounded blocks remove one NumPy allocation per subset, while
            # retaining the same inclusion-exclusion identity and O(S) storage.
            vectors=[];coefficients=[]
            def reduce_block():
                g=np.asarray(vectors);c=np.asarray(coefficients)
                covariance[:]+=g.T@(g*c[:,None])
                absolute_terms[:]+=g.T@(g*abs(c[:,None]))
            for subset,g in sums.items():
                vectors.append(g)
                coefficients.append((-1)**(len(subset)+1)*prod(1-p[v] for v in subset))
                if len(vectors)==512:
                    reduce_block();vectors.clear();coefficients.clear()
            if vectors:reduce_block()
            # BLAS summation can introduce last-bit asymmetry.
            covariance[:]=.5*covariance+.5*covariance.T
            absolute_terms[:]=.5*absolute_terms+.5*absolute_terms.T
        except FloatingPointError:raise ValueError("Covariance calculation overflows; probabilities too small") from None
    if not np.isfinite(totals).all() or not np.isfinite(covariance).all():
        raise ValueError("Nonfinite calculation")
    warnings=[]
    if np.any(counts==0):warnings.append("zero_observed_count_does_not_establish_population_absence")
    eigenvalues=np.linalg.eigvalsh(covariance)
    scale=max(1,float(abs(covariance).max()))
    if float(eigenvalues.min()) < -1e-10*scale:
        warnings.append("estimated_covariance_not_positive_semidefinite; do_not_use_for_joint_Wald_inference")
    if np.any(np.diag(covariance)<0):
        warnings.append("negative_variance_from_numerical_cancellation; inspect_before_use")
    # Finite unbiased variances do not establish asymptotic normality.
    warnings.append("normal_approximation_and_interval_coverage_not_certified")
    canonical=dict(nodes=sorted(nodes),edges=sorted(unique),probabilities=p,design=design,complete_induced=True)
    digest=hashlib.sha256(json.dumps(canonical,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return dict(software="motifsurvey-prototype",version=__version__,design=design,
                estimand="Fixed-population non-induced motif totals",motif_order=list(NAMES),
                observed_counts=counts.tolist(),estimated_totals=totals.tolist(),
                estimated_covariance=covariance.tolist(),warnings=warnings,
                diagnostics=dict(observed_nodes=len(nodes),observed_edges=len(unique),
                                 motif_instances=int(counts.sum()),subset_keys=len(sums),
                                 minimum_covariance_eigenvalue=float(eigenvalues.min()),
                                 covariance_absolute_term_sum=absolute_terms.tolist()),
                input_sha256=digest)
