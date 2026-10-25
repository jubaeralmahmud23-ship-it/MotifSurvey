"""Follow-up Monte Carlo oracle using direct weighted motif totals only.

Triggered after inspecting the original 600-round evaluation. Retain that
evaluation unchanged. This fixed 10,000-round follow-up checks true variance;
it does not provide new MotifSurvey covariance-estimation results.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "examples"))
from closure_design import load_reference, make_designs


def direct_contrast(graph, nodes, probabilities, uniforms, threshold=0.25):
    """Weighted wedge degree identity plus literal triangle vertex triples."""
    a, b, triangles = graph
    weights = (uniforms < probabilities) / probabilities
    n = len(nodes)
    neighbor = np.bincount(a, weights=weights[b], minlength=n) + np.bincount(b, weights=weights[a], minlength=n)
    neighbor_squared = np.bincount(a, weights=weights[b] ** 2, minlength=n) + np.bincount(b, weights=weights[a] ** 2, minlength=n)
    wedge = np.dot(weights, neighbor ** 2 - neighbor_squared) / 2
    triangle = np.sum(weights[triangles[:, 0]] * weights[triangles[:, 1]] * weights[triangles[:, 2]])
    return float(3 * triangle - threshold * wedge)


def prepare(graph, nodes):
    index = {v: i for i, v in enumerate(nodes)}
    edges = np.array([(index[u], index[v]) for u, v in graph.edges], dtype=int)
    adj = {v: set(graph[v]) for v in nodes}
    triangles = np.array([(index[u], index[v], index[w]) for u in nodes for v in sorted(adj[u]) if v > u for w in sorted(adj[u] & adj[v]) if w > v], dtype=int).reshape(-1, 3)
    return edges[:, 0], edges[:, 1], triangles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("data-cache"))
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    target = args.results / "variance_confirmation.json"
    if target.exists():
        raise FileExistsError(target)
    summary = json.loads((args.results / "summary.json").read_text())
    graph, labels, sources = load_reference(args.cache)
    nodes = sorted(graph)
    designs = make_designs(nodes, labels, 201)
    prepared = prepare(graph, nodes)
    # The new direct formula must first reproduce all recorded pilot/evaluation
    # contrast estimates, using their original seeds and allocations.
    raw = json.loads((args.results / "realizations.json").read_text())
    max_relative = 0.0
    for stage in raw.values():
        for name, rows in stage.items():
            for row in rows:
                uniforms = np.random.default_rng(row["seed"]).random(len(nodes))
                actual = direct_contrast(prepared, nodes, designs[name], uniforms)
                expected = row["contrast"]["estimate"]
                max_relative = max(max_relative, abs(actual - expected) / max(1, abs(expected)))
                np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-7)
    values = {name: [] for name in designs}
    for i in range(10000):
        uniforms = np.random.default_rng(2026110000 + i).random(len(nodes))
        for name, probabilities in designs.items():
            values[name].append(direct_contrast(prepared, nodes, probabilities, uniforms))
        if (i + 1) % 1000 == 0:
            print(f"Independent totals: {i + 1}/10000 paired rounds", flush=True)
    truth = summary["metadata"]["population"]["contrast"]
    result = {
        "reason": "Follow-up after initial evaluation showed Monte Carlo uncertainty and a 0.870-0.923 mean-estimated/exact variance ratio. Original pilot/evaluation retained unchanged.",
        "implementation": "Direct weighted degree identity for wedges and explicit population triangle triples; no subset covariance accumulation",
        "seeds": {"start": 2026110000, "rounds": 10000, "common_uniforms_across_designs": True},
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reproduced_original_contrasts": 2400,
        "maximum_scaled_replay_error": max_relative,
        "sources": sources, "designs": {}, "raw_contrasts": values,
        "restriction": "These are direct-total oracle calculations, not additional runs of the MotifSurvey covariance estimator or an independent human review.",
    }
    for name, estimates in values.items():
        errors = np.array(estimates) - truth
        squared = errors ** 2
        exact = summary["metadata"]["designs"][name]["exact_contrast_variance"]
        mse_mcse = squared.std(ddof=1) / np.sqrt(len(squared))
        result["designs"][name] = {
            "mean_contrast": float(np.mean(estimates)),
            "bias_in_mcse": float(errors.mean() / (errors.std(ddof=1) / np.sqrt(len(errors)))),
            "mse": float(squared.mean()), "exact_variance": exact,
            "mse_to_exact_variance": float(squared.mean() / exact),
            "mse_mcse": float(mse_mcse),
            "mse_error_in_mcse": float((squared.mean() - exact) / mse_mcse),
        }
    target.write_text(json.dumps(result, indent=2, allow_nan=False))
    print(json.dumps({k: v for k, v in result.items() if k != "raw_contrasts"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
