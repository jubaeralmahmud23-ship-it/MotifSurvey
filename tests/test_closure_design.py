"""Independent mathematical checks for the research example."""
from itertools import combinations, product
from math import prod
from pathlib import Path
import sys
import unittest

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
from closure_design import budget_probabilities, closure_contrast, exact_population_covariance
from motifsurvey import estimate
from motifsurvey import __version__
from importlib.metadata import version


def units(graph):
    result = [(frozenset(e), 0) for e in graph.edges]
    for center in graph:
        result.extend((frozenset((center, u, v)), 1) for u, v in combinations(sorted(graph[center]), 2))
    for a, b, c in combinations(sorted(graph), 3):
        if graph.has_edge(a, b) and graph.has_edge(a, c) and graph.has_edge(b, c):
            result.append((frozenset((a, b, c)), 2))
    return result


class ClosureDesignTests(unittest.TestCase):
    def test_installed_and_result_versions_agree(self):
        result = estimate([], [], {}, design="independent_node", complete_induced=True)
        self.assertEqual(result["version"], __version__)
        self.assertEqual(__version__, version("motifsurvey-prototype"))

    def test_population_oracle_against_pairs_all_four_node_graphs(self):
        possible = list(combinations(range(4), 2))
        for mask in range(64):
            graph = nx.Graph()
            graph.add_nodes_from(range(4))
            graph.add_edges_from(e for i, e in enumerate(possible) if mask >> i & 1)
            for probabilities in ({i: 0.2 for i in graph}, {i: p for i, p in enumerate((0.15, 0.4, 0.7, 1.0))}):
                reference = np.zeros((3, 3))
                rows = units(graph)
                for a, t in rows:
                    for b, u in rows:
                        reference[t, u] += 1 / prod(probabilities[v] for v in a & b) - 1
                np.testing.assert_allclose(exact_population_covariance(graph, probabilities), reference, rtol=1e-11, atol=1e-9)

    def test_contrast_expectation_variance_and_covariance_term(self):
        graph = nx.Graph([(0, 1), (1, 2), (0, 2), (2, 3)])
        probabilities = {0: 0.2, 1: 0.4, 2: 0.7, 3: 0.9}
        counts = np.bincount([t for _, t in units(graph)], minlength=3)
        coefficients = np.array([0, -0.25, 3.0])
        truth = coefficients @ counts
        mean, mse, expected_variance = 0.0, 0.0, 0.0
        for bits in product((0, 1), repeat=4):
            sample = graph.subgraph([v for v in graph if bits[v]])
            probability = prod(probabilities[v] if bits[v] else 1 - probabilities[v] for v in graph)
            result = estimate([str(v) for v in sample], [(str(u), str(v)) for u, v in sample.edges], {str(v): probabilities[v] for v in sample}, design="independent_node", complete_induced=True)
            contrast = closure_contrast(result)
            mean += probability * contrast["estimate"]
            mse += probability * (contrast["estimate"] - truth) ** 2
            expected_variance += probability * contrast["variance"]
        exact = coefficients @ exact_population_covariance(graph, probabilities) @ coefficients
        self.assertAlmostEqual(mean, truth, places=10)
        self.assertAlmostEqual(mse, exact, places=9)
        self.assertAlmostEqual(expected_variance, exact, places=9)

    def test_probability_budget(self):
        values = budget_probabilities([1, 2, 4, 20], 1.2)
        self.assertAlmostEqual(sum(values), 1.2, places=12)
        self.assertTrue(np.all((values >= 0.05) & (values <= 0.8)))
        with self.assertRaises(ValueError):
            budget_probabilities([0, 1], 0.4)

    def test_negative_contrast_variance_preserved(self):
        result = {"estimated_totals": [0, 4, 1], "estimated_covariance": [[0, 0, 0], [0, 1, 100], [0, 100, 1]]}
        self.assertLess(closure_contrast(result)["variance"], 0)


if __name__ == "__main__":
    unittest.main()
