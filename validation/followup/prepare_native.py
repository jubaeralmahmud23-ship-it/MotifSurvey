import json
from pathlib import Path
import networkx as nx
import numpy as np
from stress import calc
B=Path(__file__).resolve().parent
fixtures=[];expected=[]
for case in range(24):
 g=nx.fast_gnp_random_graph(6+case%5,.15+.12*(case%5),seed=2000+case)
 if g.number_of_edges()<2:g.add_edges_from([(0,1),(1,2)])
 p={v:([.1,.3,.8,1.][(v+case)%4] if case%3 else .4) for v in g}
 fixtures.append(dict(name=f'correctness_{case}',n=len(g),edges=list(g.edges),p=list(p.values()),methods=['sampling','survey'],repeats=1))
 expected.append(dict(name=f'correctness_{case}',result=calc(g,p)))
for name,g in [('ring_500',nx.cycle_graph(500)),('hub_50',nx.star_graph(49)),('random_200',nx.fast_gnp_random_graph(200,.008,seed=44)),('clique_12',nx.complete_graph(12)),('ring_2000',nx.cycle_graph(2000)),('clique_25',nx.complete_graph(25))]:
 p={v:.3+.6*(v%7)/6 for v in g}
 fixtures.append(dict(name=name,n=len(g),edges=list(g.edges),p=list(p.values()),methods=['survey'] if name in ('ring_2000','clique_25') else ['sampling','survey'],repeats=3))
 expected.append(dict(name=name,result=calc(g,p)))
(B/'native_fixtures.json').write_text(json.dumps(fixtures))
(B/'native_expected.json').write_text(json.dumps(expected))
print('Prepared 24 correctness and 6 benchmark fixtures')
