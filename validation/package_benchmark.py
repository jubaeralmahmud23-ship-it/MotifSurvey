"""End-to-end installed package vs local formula references; no native R claim."""
import json,time,tracemalloc
from pathlib import Path
from itertools import combinations
import numpy as np
from motifsurvey import estimate
from preflight import instances,indexed,vectorized
B=Path(__file__).resolve().parent
rng=np.random.default_rng(440)
fixtures=[('ring',2000,[(i,(i+1)%2000) for i in range(2000)]),('hub_all_motifs',50,[(0,i) for i in range(1,50)]),('random',700,[e for e in combinations(range(700),2) if rng.random()<.004])]
results=[]
for name,n,edges in fixtures:
 p={i:float(.3+.6*(i%7)/6) for i in range(n)}
 def package():
  r=estimate([str(v) for v in range(n)],[(str(u),str(v)) for u,v in edges],{str(v):p[v] for v in range(n)},design='independent_node',complete_induced=True)
  return np.array(r['estimated_totals']),np.array(r['estimated_covariance'])
 def baseline(fn):return lambda:fn(instances(list(range(n)),edges),p)
 record=dict(name=name,nodes=n,edges=len(edges),motif_instances=len(instances(list(range(n)),edges)),timing_scope='End-to-end, includes motif enumeration and package validation/hash; references omit validation/hash')
 for label,fn in [('installed_package',package),('indexed_reference',baseline(indexed)),('dense_reference',baseline(vectorized))]:
  timings=[]
  for _ in range(3):
   t=time.perf_counter();a,c=fn();timings.append(time.perf_counter()-t)
  tracemalloc.start();fn();_,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
  record[label]=dict(seconds=timings,median_seconds=float(np.median(timings)),peak_traced_bytes=peak,totals=a.tolist(),covariance=c.tolist())
 for label in ('indexed_reference','dense_reference'):
  np.testing.assert_allclose(record[label]['covariance'],record['installed_package']['covariance'],rtol=1e-8,atol=1e-8)
 results.append(record);(B/'package_benchmark.json').write_text(json.dumps(results,indent=2));print(json.dumps(record),flush=True)
