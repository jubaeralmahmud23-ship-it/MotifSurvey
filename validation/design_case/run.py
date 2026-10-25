"""Execute the frozen closure-design protocol and retain every realization."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time

import networkx as nx
import numpy as np
import motifsurvey
import motifsurvey.core

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "examples"))
from closure_design import load_reference, make_designs, sampled_estimate, closure_contrast, exact_population_covariance


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")


def run(cache, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError("Choose a new results directory to preserve previous evidence")
    output.mkdir(parents=True)
    protocol_path = Path(__file__).with_name("protocol.json")
    protocol = json.loads(protocol_path.read_text())
    graph, labels, sources = load_reference(cache)
    nodes = sorted(graph)
    designs = make_designs(nodes, labels, protocol["expected_nodes"])
    truth = np.array([graph.number_of_edges(), sum(d * (d - 1) // 2 for _, d in graph.degree()), sum(nx.triangles(graph).values()) // 3], dtype=float)
    coefficients = np.array([0, -protocol["threshold"], 3.0])
    target = float(coefficients @ truth)
    metadata = {
        "protocol": protocol,
        "protocol_sha256": hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
        "sources": sources,
        "population": {"nodes": len(graph), "edges": graph.number_of_edges(), "departments": len(set(labels.values())), "totals": truth.tolist(), "transitivity": float(3 * truth[2] / truth[1]), "contrast": target},
        "environment": {"python": platform.python_version(), "platform": platform.platform(), "numpy": np.__version__, "networkx": nx.__version__, "motifsurvey": motifsurvey.__version__},
        "core_sha256": hashlib.sha256(Path(motifsurvey.core.__file__).read_bytes()).hexdigest(),
        "designs": {},
    }
    for name, probabilities in designs.items():
        pmap = dict(zip(nodes, probabilities))
        covariance = exact_population_covariance(graph, pmap)
        metadata["designs"][name] = {
            "probabilities": probabilities.tolist(), "node_order": nodes,
            "expected_nodes": float(sum(probabilities)),
            "expected_edges": float(sum(pmap[u] * pmap[v] for u, v in graph.edges)),
            "probability_range": [float(probabilities.min()), float(probabilities.max())],
            "exact_covariance": covariance.tolist(),
            "exact_contrast_variance": float(coefficients @ covariance @ coefficients),
            "exact_diagonal_only_variance": float((coefficients ** 2) @ np.diag(covariance)),
        }
    write(output / "metadata.json", metadata)
    raw = {stage: {name: [] for name in designs} for stage in ("pilot", "evaluation")}
    selections = None
    for stage in ("pilot", "evaluation"):
        count = protocol[f"{stage}_repetitions_per_design"]
        start_seed = protocol[f"{stage}_seed_start"]
        for i in range(count):
            seed = start_seed + i
            uniforms = np.random.default_rng(seed).random(len(nodes))
            for name, probabilities in designs.items():
                begin = time.perf_counter()
                result = sampled_estimate(graph, nodes, probabilities, uniforms)
                contrast = closure_contrast(result, protocol["threshold"])
                raw[stage][name].append({
                    "seed": seed, "selected_nodes": int((uniforms < probabilities).sum()),
                    "observed_counts": result["observed_counts"],
                    "estimated_totals": result["estimated_totals"],
                    "estimated_covariance": result["estimated_covariance"],
                    "contrast": contrast, "warnings": result["warnings"],
                    "seconds": time.perf_counter() - begin,
                })
            if (i + 1) % 50 == 0:
                write(output / "realizations.json", raw)
                print(f"{stage}: {i + 1}/{count} paired rounds ({len(designs)} designs)", flush=True)
        if stage == "pilot":
            full = {name: float(np.mean([r["contrast"]["variance"] for r in rows])) for name, rows in raw[stage].items()}
            diagonal = {name: float(np.mean([r["contrast"]["variance_ignoring_covariance"] for r in rows])) for name, rows in raw[stage].items()}
            selections = {"mean_estimated_variances": full, "mean_diagonal_only_variances": diagonal, "full_covariance_choice": min(full, key=full.get), "diagonal_only_choice": min(diagonal, key=diagonal.get)}
            write(output / "pilot_selection.json", selections)
            print(json.dumps(selections), flush=True)
    summary = {"metadata": metadata, "pilot_selection": selections, "evaluation": {}}
    uniform_errors = np.array([r["contrast"]["estimate"] - target for r in raw["evaluation"]["uniform"]])
    bootstrap = np.random.default_rng(2026104500).integers(0, len(uniform_errors), size=(2000, len(uniform_errors)))
    for name, rows in raw["evaluation"].items():
        errors = np.array([r["contrast"]["estimate"] - target for r in rows])
        variance_estimates = np.array([r["contrast"]["variance"] for r in rows])
        empirical_variance = float(np.var(errors, ddof=1))
        mean_estimated_variance = float(variance_estimates.mean())
        ratio_boot = np.mean(errors[bootstrap] ** 2, axis=1) / np.mean(uniform_errors[bootstrap] ** 2, axis=1)
        exact = metadata["designs"][name]["exact_contrast_variance"]
        diagonal = np.array([r["contrast"]["variance_ignoring_covariance"] for r in rows])
        summary["evaluation"][name] = {
            "repetitions": len(rows), "mean_contrast": float(target + errors.mean()),
            "bias_in_mcse": float(errors.mean() / np.sqrt(empirical_variance / len(rows))),
            "empirical_variance": empirical_variance, "empirical_mse": float(np.mean(errors ** 2)),
            "exact_variance": exact, "mean_estimated_variance": mean_estimated_variance,
            "mean_estimated_to_exact_variance": mean_estimated_variance / exact,
            "estimated_variance_mean_mcse": float(variance_estimates.std(ddof=1) / np.sqrt(len(rows))),
            "estimated_variance_mean_error_in_mcse": float((mean_estimated_variance - exact) / (variance_estimates.std(ddof=1) / np.sqrt(len(rows)))),
            "mean_estimated_to_empirical_variance": mean_estimated_variance / empirical_variance,
            "exact_variance_ratio_to_uniform": exact / metadata["designs"]["uniform"]["exact_contrast_variance"],
            "mse_ratio_to_uniform": float(np.mean(errors ** 2) / np.mean(uniform_errors ** 2)),
            "mse_ratio_paired_bootstrap_95": np.quantile(ratio_boot, [0.025, 0.975]).tolist(),
            "mean_diagonal_only_to_exact_variance": float(diagonal.mean() / exact),
            "mean_nodes": float(np.mean([r["selected_nodes"] for r in rows])),
            "mean_edges": float(np.mean([r["observed_counts"][0] for r in rows])),
            "negative_estimated_variance": int((variance_estimates < 0).sum()),
            "indefinite_covariance": sum(any("not_positive_semidefinite" in warning for warning in r["warnings"]) for r in rows),
        }
    summary["interpretation"] = "Bootstrap intervals quantify simulation uncertainty in the paired MSE ratio, not confidence intervals for the population contrast. Pilot and evaluation seeds are disjoint; designs share uniforms at each seed. Exact population values are validation oracles, not sample-accessible inputs."
    write(output / "realizations.json", raw)
    write(output / "summary.json", summary)
    print(json.dumps({"pilot_selection": selections, "evaluation": summary["evaluation"]}, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("data-cache"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.cache, args.output)
