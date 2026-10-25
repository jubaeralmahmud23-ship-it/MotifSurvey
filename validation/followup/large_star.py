"""Higher explicit resource budget; exact closed-form oracle for a star."""
import json,math,subprocess,sys,time
from pathlib import Path
import numpy as np
import networkx as nx
from benchmark import peak_memory
from stress import calc
B=Path(__file__).resolve().parent
L=999;g=nx.star_graph(L);p=.5;W=math.comb(L,2)
start=time.perf_counter();r=calc(g,dict.fromkeys(g,p),max_instances=500_000,max_subsets=1_000_000);seconds=time.perf_counter()-start
vE=(L*(1-p*p)+L*(L-1)*(1-p))/p**4
vW=(W*(1-p**3)+W*2*(L-2)*(1-p*p)+W*math.comb(L-2,2)*(1-p))/p**6
cEW=(2*W*(1-p*p)+(L-2)*W*(1-p))/p**5
expected=np.array([[vE,cEW,0],[cEW,vW,0],[0,0,0]])
np.testing.assert_allclose(r['estimated_covariance'],expected,rtol=1e-12,atol=1e-8)
np.testing.assert_allclose(r['estimated_totals'],[L/p**2,W/p**3,0],rtol=1e-12,atol=1e-8)
out=dict(version=r['version'],nodes=L+1,instances=L+W,subset_keys=r['diagnostics']['subset_keys'],explicit_max_instances=500000,explicit_max_subsets=1000000,seconds=seconds,process_peak_working_set=peak_memory(),closed_form_covariance=expected.tolist(),actual=r,all_passed=True,comparison='No quadratic baseline executed at this scale. Default instance budget correctly rejects this input; this run deliberately raised it.')
(B/'large_star_results.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='actual'}),flush=True)
