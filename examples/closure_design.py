"""Compare independent-vertex sampling designs for a closure contrast.

This is a controlled design-planning experiment on a reference population,
not an assertion about the mechanism that originally collected email data.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import gzip
import hashlib
import json
import urllib.request

import networkx as nx
import numpy as np

from motifsurvey import estimate

GRAPH_URL = "https://snap.stanford.edu/data/email-Eu-core.txt.gz"
GRAPH_SHA256 = "4b47acdb80197b085fe63c819c357ae488131ee904ed93d1b219a68b0f9e245f"
LABEL_URL = "https://snap.stanford.edu/data/email-Eu-core-department-labels.txt.gz"
LABEL_SHA256 = "e5abe5b4581a480032a63adcf2576c161785f45692642c6ebb0b1276f0f33669"


def fetch_verified(url, target, expected_sha256=None):
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        request = urllib.request.Request(url, headers={"User-Agent": "MotifSurvey research example"})
        with urllib.request.urlopen(request, timeout=60) as response:
            target.write_bytes(response.read())
    data = target.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError(f"Dataset checksum mismatch: {target.name}")
    return data, digest


def load_reference(cache, label_sha256=LABEL_SHA256):
    raw, graph_hash = fetch_verified(GRAPH_URL, Path(cache) / "email-Eu-core.txt.gz", GRAPH_SHA256)
    labels_raw, labels_hash = fetch_verified(LABEL_URL, Path(cache) / "email-Eu-core-department-labels.txt.gz", label_sha256)
    graph = nx.Graph()
    for line in gzip.decompress(raw).decode().splitlines():
        if line and not line.startswith("#"):
            u, v = map(int, line.split())
            graph.add_nodes_from((u, v))
            if u != v:
                graph.add_edge(u, v)
    labels = dict(tuple(map(int, line.split())) for line in gzip.decompress(labels_raw).decode().splitlines() if line and not line.startswith("#"))
    if set(labels) != set(graph):
        raise ValueError("Department roster must cover exactly the reference nodes")
    return graph, labels, {"graph_sha256": graph_hash, "labels_sha256": labels_hash}


def budget_probabilities(scores, expected_nodes, lower=0.05, upper=0.8):
    """Calibrate bounded independent selection probabilities to an expected cost."""
    scores = np.asarray(scores, dtype=float)
    if not np.all(np.isfinite(scores)) or np.any(scores <= 0):
        raise ValueError("Allocation scores must be finite and positive")
    if not len(scores) * lower <= expected_nodes <= len(scores) * upper:
        raise ValueError("Expected sample size is outside the probability bounds")
    lo, hi = 0.0, upper / scores.min()
    for _ in range(100):
        mid = (lo + hi) / 2
        if np.clip(mid * scores, lower, upper).sum() < expected_nodes:
            lo = mid
        else:
            hi = mid
    return np.clip((lo + hi) * scores / 2, lower, upper)


def make_designs(nodes, labels, expected_nodes):
    sizes = Counter(labels.values())
    group_sizes = np.array([sizes[labels[v]] for v in nodes], dtype=float)
    return {
        "uniform": np.full(len(nodes), expected_nodes / len(nodes)),
        "small_department": budget_probabilities(1 / np.sqrt(group_sizes), expected_nodes),
        "large_department": budget_probabilities(np.sqrt(group_sizes), expected_nodes),
    }


def sampled_estimate(graph, nodes, probabilities, uniforms):
    selected = [v for v, p, u in zip(nodes, probabilities, uniforms) if u < p]
    pmap = dict(zip(nodes, probabilities))
    sample = graph.subgraph(selected)
    result = estimate(
        [str(v) for v in selected],
        [(str(u), str(v)) for u, v in sample.edges],
        {str(v): float(pmap[v]) for v in selected},
        design="independent_node", complete_induced=True,
    )
    return result


def closure_contrast(result, threshold=0.25):
    """Unbiased 3*T_triangle - threshold*T_wedge and estimated variance.

    Its sign is the sign of transitivity minus the threshold when the
    population wedge count is positive. It is not a ratio estimator or a test.
    An individual estimated variance can be negative and is not clipped.
    """
    coefficients = np.array([0.0, -threshold, 3.0])
    totals = np.asarray(result["estimated_totals"])
    covariance = np.asarray(result["estimated_covariance"])
    return {
        "estimate": float(coefficients @ totals),
        "variance": float(coefficients @ covariance @ coefficients),
        "variance_ignoring_covariance": float((coefficients ** 2) @ np.diag(covariance)),
    }


def exact_population_covariance(graph, probabilities):
    """Independent full-population oracle; not available from a sample.

    Covariance = sum_A product(1/p_v-1) h_A h_A', where h_A counts
    population motifs containing A. Specialized incidence counts avoid
    storing all triple keys. Tests compare this with literal ordered pairs
    and exhaustive sample distributions.
    """
    nodes = sorted(graph)
    adj = {v: set(graph[v]) for v in nodes}
    degrees = dict(graph.degree())
    triangle_counts = nx.triangles(graph)
    q = {v: 1 / probabilities[v] - 1 for v in nodes}
    covariance = np.zeros((3, 3))
    for v in nodes:
        d = degrees[v]
        incidence = np.array([d, d * (d - 1) // 2 + sum(degrees[u] - 1 for u in adj[v]), triangle_counts[v]], dtype=float)
        covariance += q[v] * np.outer(incidence, incidence)
    for u, v in combinations(nodes, 2):
        edge = int(v in adj[u])
        common = len(adj[u] & adj[v])
        wedge = common + edge * (degrees[u] + degrees[v] - 2)
        if edge or wedge:
            incidence = np.array([edge, wedge, edge * common], dtype=float)
            covariance += q[u] * q[v] * np.outer(incidence, incidence)
    open_triples = 0.0
    triangle_triples = 0.0
    for center in nodes:
        for u, v in combinations(sorted(adj[center]), 2):
            value = q[center] * q[u] * q[v]
            if v not in adj[u]:
                open_triples += value
            elif center < u < v:
                triangle_triples += value
    covariance[1, 1] += open_triples + 9 * triangle_triples
    covariance[1, 2] += 3 * triangle_triples
    covariance[2, 1] += 3 * triangle_triples
    covariance[2, 2] += triangle_triples
    return covariance


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON MotifSurvey result")
    parser.add_argument("--threshold", type=float, default=0.25)
    args = parser.parse_args()
    print(json.dumps(closure_contrast(json.loads(args.input.read_text()), args.threshold), indent=2))
