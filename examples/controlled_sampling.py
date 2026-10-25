"""A controlled experiment. This does not model unknown real-data missingness.

Install the validation extra for NetworkX. No file is written by this example.
"""
import json
import networkx as nx
import numpy as np
from motifsurvey import estimate

population = nx.barabasi_albert_graph(800, 3, seed=190)
# We control this artificial observation design, so probabilities are known.
p = {v: .15 + .3 * (v % 7) / 6 for v in population}
rng = np.random.default_rng(191)
observed = population.subgraph([v for v in population if rng.random() < p[v]])
result = estimate(
    nodes=[str(v) for v in observed],
    edges=[(str(u), str(v)) for u, v in observed.edges],
    probabilities={str(v): p[v] for v in observed},
    design="independent_node",
    complete_induced=True,
)
truth = [population.number_of_edges(),
         sum(d * (d - 1) // 2 for _, d in population.degree()),
         sum(nx.triangles(population).values()) // 3]
print(json.dumps({"population_truth_for_this_simulation_only": truth,
                  "sample_analysis": result}, indent=2, allow_nan=False))
