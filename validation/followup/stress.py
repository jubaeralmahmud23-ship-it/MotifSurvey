"""Independent Decimal pair oracle against the installed MotifSurvey."""
import json,sys,time,platform
from pathlib import Path
from decimal import Decimal,localcontext
from itertools import combinations
import numpy as np
import networkx as nx
import motifsurvey
from motifsurvey import estimate,ResourceLimitError
B=Path(__file__).resolve().parent
sys.path.insert(0,str(B.parent/'motif_search' if (B.parent/'motif_search').exists() else B.parent))
from preflight import instances

def calc(g,p,**kwargs):
 return estimate([str(v) for v in g],[(str(u),str(v)) for u,v in g.edges],{str(v):p[v] for v in g},design='independent_node',complete_induced=True,**kwargs)

def decimal_oracle(g,p):
 rows=instances(list(g),list(g.edges))
 with localcontext() as ctx:
  ctx.prec=90
  pd={v:Decimal.from_float(float(x)) for v,x in p.items()}
  def product(s):
   value=Decimal(1)
   for v in s:value*=pd[v]
   return value
  weights=[1/product(s) for s,t in rows];cov=[[Decimal(0) for _ in range(3)] for _ in range(3)]
  for i,(s,t) in enumerate(rows):
   for j,(a,b) in enumerate(rows):cov[t][b]+=weights[i]*weights[j]*(1-product(s&a))
  return np.array([[float(x) for x in r] for r in cov])

if __name__=='__main__':
 graphs={'clique':nx.complete_graph(6),'star':nx.star_graph(8),'cycle':nx.cycle_graph(8),'bipartite':nx.complete_bipartite_graph(4,4),'disconnected':nx.disjoint_union(nx.complete_graph(4),nx.path_graph(4))}
 regimes={
  'census':lambda v:1.,
  'almost_census':lambda v:float(np.nextafter(1.,0.)),
  'uniform_low':lambda v:1e-20,
  'heterogeneous_extreme':lambda v:[1e-12,.02,.8,1.][v%4],
  'ordinary':lambda v:.15+.8*(v%7)/6,
 }
 results=[];failures=[]
 for name,g in graphs.items():
  for regime,pfn in regimes.items():
   p={v:pfn(v) for v in g};r=calc(g,p);actual=np.array(r['estimated_covariance']);oracle=decimal_oracle(g,p)
   error=float(np.max(abs(actual-oracle)/np.maximum(abs(oracle),1e-300)))
   zero=float(np.max(abs(actual[oracle==0]),initial=0))
   good=bool(error<=1e-8 and zero==0)
   item=dict(graph=name,regime=regime,maximum_relative_error=error,maximum_error_on_exact_zeros=zero,passed=good,warnings=r['warnings']);results.append(item)
   if not good:failures.append(item)
   print(json.dumps(item),flush=True)
 limits=[]
 for name,g,pv,kw in [('hub_1000_default',nx.star_graph(999),.5,{}),('clique_100_default',nx.complete_graph(100),.5,{}),('inverse_overflow',nx.complete_graph(3),1e-200,{}),('covariance_overflow',nx.complete_graph(3),1e-60,{}),('subset_cap',nx.cycle_graph(20),.5,dict(max_subsets=5))]:
  start=time.perf_counter()
  try:calc(g,{v:pv for v in g},**kw);result='unexpected success'
  except (ValueError,ResourceLimitError) as e:result=type(e).__name__+': '+str(e)
  limits.append(dict(case=name,result=result,seconds=time.perf_counter()-start))
 out=dict(oracle='90-digit Decimal direct ordered motif pairs; same IEEE input probabilities',cases=results,failures=failures,resource_checks=limits,versions=dict(motifsurvey=motifsurvey.__version__,python=sys.version,numpy=np.__version__,networkx=nx.__version__,platform=platform.platform()))
 (B/'stress_results.json').write_text(json.dumps(out,indent=2,allow_nan=False))
 print(json.dumps(dict(cases=len(results),failures=len(failures),resource_checks=limits)),flush=True)
 if failures or any(x['result']=='unexpected success' for x in limits):
  raise SystemExit('High-precision or expected-rejection validation failed')
