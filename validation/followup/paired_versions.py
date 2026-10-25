import json,time
from pathlib import Path
import numpy as np
from motifsurvey import estimate
from core_v010 import estimate as old_estimate
from benchmark import fixture
B=Path(__file__).resolve().parent
results=[];rng=np.random.default_rng(555)
for name in ['ring_2000','hub_50','random_700','clique_25','hub_200','heavy_tail_1500','ring_10000']:
 g=fixture(name);payload=dict(nodes=[str(v) for v in g],edges=[(str(u),str(v)) for u,v in g.edges],probabilities={str(v):.3+.6*(v%7)/6 for v in g},design='independent_node',complete_induced=True,max_instances=2000000,max_subsets=4000000)
 funcs={'0.1.0':old_estimate,'0.2.0':estimate};times={v:[] for v in funcs}
 for f in funcs.values():f(**payload)
 outputs={}
 for repeat in range(7):
  for version in rng.permutation(list(funcs)):
   start=time.perf_counter();outputs[version]=funcs[version](**payload);times[version].append(time.perf_counter()-start)
 np.testing.assert_allclose(outputs['0.1.0']['estimated_covariance'],outputs['0.2.0']['estimated_covariance'],rtol=1e-8,atol=1e-8)
 row=dict(fixture=name,seconds=times,old_median=float(np.median(times['0.1.0'])),new_median=float(np.median(times['0.2.0'])))
 row['old_over_new']=row['old_median']/row['new_median'];results.append(row)
 (B/'paired_version_results.json').write_text(json.dumps(dict(scope='Seven repetitions, random version order within each repetition, one warm-up each. Same prepared Python inputs, both include validation/hash. Same process per all fixtures; no OS memory comparison in this experiment.',records=results),indent=2))
 print(json.dumps(row),flush=True)
