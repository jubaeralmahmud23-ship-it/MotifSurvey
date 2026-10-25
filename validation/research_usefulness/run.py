"""Reference-budget sensitivity, archived-sample reanalysis, and R fixtures.

MIT License. Requires the released MotifSurvey 0.2.1 and its validation extra.
"""
import argparse
from itertools import combinations
import hashlib
import json
from pathlib import Path
import platform
import sys

import networkx as nx
import numpy as np
import motifsurvey
import motifsurvey.core

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "examples"))
from closure_design import load_reference, make_designs, sampled_estimate, exact_population_covariance
from closure_profile import closure_profile


def expected_observed_counts(graph, p):
    """Exact first moments, using weighted degrees and triangle triples."""
    edges = sum(p[u] * p[v] for u, v in graph.edges)
    wedges = sum(p[v] * ((sum(p[u] for u in graph[v])) ** 2 - sum(p[u] ** 2 for u in graph[v])) / 2 for v in graph)
    triangles = sum(p[u] * p[v] * p[w] for u in graph for v, w in combinations(sorted(graph[u]), 2) if u < v < w and graph.has_edge(v, w))
    return np.array([edges, wedges, triangles])


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def run(cache, output):
    output.mkdir(parents=True, exist_ok=False)
    protocol_file = Path(__file__).with_name("protocol.json")
    protocol = json.loads(protocol_file.read_text())
    graph, labels, sources = load_reference(cache)
    nodes = sorted(graph)
    truth = np.array([graph.number_of_edges(), sum(d * (d - 1) // 2 for _, d in graph.degree()), sum(nx.triangles(graph).values()) // 3], dtype=float)
    original_path = ROOT / "validation/design_case/results/realizations.json"
    original = json.loads(original_path.read_text())
    original_meta = json.loads((original_path.parent / "metadata.json").read_text())
    rows, reanalysis, fixtures, expected = [], {}, [], []
    for budget in protocol["expected_vertex_budgets"]:
        designs = make_designs(nodes, labels, budget)
        pbar = budget / len(nodes)
        for name, p in designs.items():
            pmap = dict(zip(nodes, p))
            covariance = exact_population_covariance(graph, pmap)
            observed_mean = expected_observed_counts(graph, pmap)
            uniform_correction_mean = observed_mean / pbar ** np.array([2, 3, 3])
            if name == "uniform":
                np.testing.assert_allclose(uniform_correction_mean, truth, rtol=1e-10, atol=1e-8)
            if budget == 201:
                np.testing.assert_allclose(covariance, original_meta["designs"][name]["exact_covariance"], rtol=1e-12, atol=1e-6)
                saved = original["evaluation"][name]
                methods = {"ht": np.array([x["estimated_totals"] for x in saved]),
                           "mean_probability_correction": np.array([x["observed_counts"] for x in saved]) / pbar ** np.array([2, 3, 3])}
                reanalysis[name] = {"realizations": len(saved), "methods": {}}
                for method, totals in methods.items():
                    mean = totals.mean(axis=0)
                    mcse = totals.std(axis=0, ddof=1) / np.sqrt(len(totals))
                    c = np.array([0, -.25, 3])
                    contrasts = totals @ c
                    target = truth @ c
                    reanalysis[name]["methods"][method] = {
                        "mean_totals": mean.tolist(), "relative_mean_error_percent": (100 * (mean / truth - 1)).tolist(),
                        "mean_mcse": mcse.tolist(), "mean_contrast": float(contrasts.mean()),
                        "contrast_mean_error_mcse": float((contrasts.mean() - target) / (contrasts.std(ddof=1) / np.sqrt(len(totals)))),
                        "contrast_mse": float(np.mean((contrasts - target) ** 2)),
                    }
                uniforms = np.random.default_rng(protocol["native_comparison_seed"]).random(len(nodes))
                result = sampled_estimate(graph, nodes, p, uniforms)
                selected = [v for v, prob, u in zip(nodes, p, uniforms) if u < prob]
                mapping = {v: i for i, v in enumerate(selected)}
                fixtures.append({"name": name, "n": len(selected), "edges": [[mapping[u], mapping[v]] for u, v in graph.subgraph(selected).edges], "p": [pmap[v] for v in selected], "methods": ["survey"], "repeats": 1})
                expected.append({"name": name, "result": result, "profile": closure_profile(result, protocol["closure_thresholds"])})
                np.testing.assert_allclose(result["estimated_totals"], saved[0]["estimated_totals"], rtol=1e-12)
            profiles = []
            for threshold in protocol["closure_thresholds"]:
                c = np.array([0, -threshold, 3])
                target = float(c @ truth)
                variance = float(c @ covariance @ c)
                profiles.append({"threshold": threshold, "target": target, "exact_variance": variance,
                                 "exact_standard_deviation": float(np.sqrt(variance)),
                                 "exact_diagonal_only_variance": float((c * c) @ np.diag(covariance)),
                                 "mean_probability_correction_expectation": float(c @ uniform_correction_mean)})
            rows.append({"expected_nodes": budget, "allocation": name, "probability_range": [float(p.min()), float(p.max())],
                         "expected_observed_counts": observed_mean.tolist(), "exact_covariance": covariance.tolist(),
                         "mean_probability_correction_relative_bias_percent": (100 * (uniform_correction_mean / truth - 1)).tolist(),
                         "profiles": profiles})
            print(f"Completed budget {budget}, {name}", flush=True)
    out = {"protocol": protocol, "protocol_sha256": hashlib.sha256(protocol_file.read_bytes()).hexdigest(),
           "sources": sources, "original_realizations_sha256": hashlib.sha256(original_path.read_bytes()).hexdigest(),
           "environment": {"python": platform.python_version(), "numpy": np.__version__, "networkx": nx.__version__, "motifsurvey": motifsurvey.__version__, "platform": platform.platform()},
           "core_sha256": hashlib.sha256(Path(motifsurvey.core.__file__).read_bytes()).hexdigest(),
           "population_totals": truth.tolist(), "sensitivity": rows, "archived_reanalysis": reanalysis}
    dump(output / "summary.json", out)
    dump(output / "native_fixtures.json", fixtures)
    dump(output / "native_expected.json", expected)
    print("Wrote sensitivity, all archived reanalyses, and three fixed-seed native fixtures", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("data-cache"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.cache, args.output)
