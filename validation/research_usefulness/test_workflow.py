"""Independent checks of the research companion; MIT License."""
import itertools
import json
from pathlib import Path
import sys
import unittest

import networkx as nx
import numpy as np
from motifsurvey import estimate

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "examples"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from closure_profile import closure_profile
from run import expected_observed_counts


class ResearchWorkflowTests(unittest.TestCase):
    def test_expected_counts_against_all_samples(self):
        # Directly average NetworkX observed counts over each entire sampling law.
        # All 64 four-vertex graphs, with two unequal-probability regimes.
        pairs = list(itertools.combinations(range(4), 2))
        for mask in range(64):
            graph = nx.Graph()
            graph.add_nodes_from(range(4))
            # Reverse insertion order guards against unsorted-neighbor triangle loss.
            graph.add_edges_from(edge for i, edge in reversed(list(enumerate(pairs))) if mask & (1 << i))
            for probs in ([.1, .3, .7, 1.], [.05, .8, .25, .6]):
                mean = np.zeros(3)
                for selected in itertools.product([False, True], repeat=4):
                    probability = np.prod([probs[v] if selected[v] else 1-probs[v] for v in range(4)])
                    sample = graph.subgraph([v for v in range(4) if selected[v]])
                    counts = np.array([sample.number_of_edges(), sum(d*(d-1)//2 for _, d in sample.degree()), sum(nx.triangles(sample).values())/3])
                    mean += probability * counts
                np.testing.assert_allclose(expected_observed_counts(graph, dict(enumerate(probs))), mean, rtol=1e-12, atol=1e-12)

    def test_profile_retains_negative_variance_and_provenance(self):
        result = {"motif_order": ["edges", "wedges", "triangles"], "estimated_totals": [2, 4, 1], "estimated_covariance": [[0, 0, 0], [0, 1, 2], [0, 2, 1]], "input_sha256": "test-hash", "warnings": ["indefinite"], "version": "0.2.1"}
        profile = closure_profile(result, [0, .5, 1])
        self.assertEqual(profile["input_sha256"], "test-hash")
        self.assertEqual(profile["source_warnings"], ["indefinite"])
        self.assertEqual(profile["profile"][-1]["estimate"], -1)
        self.assertEqual(profile["profile"][-1]["estimated_variance"], -2)
        self.assertTrue(profile["profile"][-1]["negative_estimated_variance"])

    def test_profile_rejects_wrong_order_and_invalid_thresholds(self):
        result = estimate(["a", "b"], [("a", "b")], {"a": .5, "b": .5}, design="independent_node", complete_induced=True)
        for thresholds in ([], [-.1], [1.1], [float("nan")], [True], [".25"], [np.bool_(True)], [np.complex128(.25+.2j)], [b'.25'], [10**5000]):
            with self.assertRaises(ValueError):
                closure_profile(result, thresholds)
        result["motif_order"] = ["triangle", "wedge", "edge"]
        with self.assertRaises(ValueError):
            closure_profile(result, [.25])

    def test_archived_native_workflow(self):
        base = Path(__file__).with_name("results")
        actual = json.loads((base / "native_results.json").read_text())
        expected = {r["name"]: r for r in json.loads((base / "native_expected.json").read_text())}
        self.assertEqual(len(actual["results"]), 3)
        for native in actual["results"]:
            target = expected[native["name"]]
            np.testing.assert_allclose(native["totals"], target["result"]["estimated_totals"], rtol=1e-8, atol=1e-8)
            np.testing.assert_allclose(native["covariance"], target["result"]["estimated_covariance"], rtol=1e-8, atol=1e-8)
            result = dict(target["result"], estimated_totals=native["totals"], estimated_covariance=native["covariance"])
            profile = closure_profile(result, [.15, .25, .35])
            for x, y in zip(profile["profile"], target["profile"]["profile"]):
                np.testing.assert_allclose([x["estimate"], x["estimated_variance"]], [y["estimate"], y["estimated_variance"]], rtol=1e-8, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
