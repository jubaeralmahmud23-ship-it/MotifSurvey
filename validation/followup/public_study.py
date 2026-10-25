"""Synthetic known-design thinning, not inference about original collection."""
import gzip,hashlib,json,math,time,urllib.request
from pathlib import Path
import networkx as nx
import numpy as np
import motifsurvey
from stress import calc
B=Path(__file__).resolve().parent
url='https://snap.stanford.edu/data/email-Eu-core.txt.gz'
path=B/'email-Eu-core.txt.gz'
if not path.exists():
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'SoftwareX-research-validation'}),timeout=30) as r:path.write_bytes(r.read())
raw=path.read_bytes()
if hashlib.sha256(raw).hexdigest()!='4b47acdb80197b085fe63c819c357ae488131ee904ed93d1b219a68b0f9e245f':
 raise ValueError('Public dataset checksum mismatch; inspect the source before running')
g=nx.Graph();records=0;loops=0
for line in gzip.decompress(raw).decode().splitlines():
 if not line or line.startswith('#'):continue
 u,v=map(int,line.split());records+=1;g.add_nodes_from((u,v))
 if u==v:loops+=1;continue
 g.add_edge(u,v)
target=np.array([g.number_of_edges(),sum(d*(d-1)//2 for _,d in g.degree()),sum(nx.triangles(g).values())//3],float)
out=dict(version=motifsurvey.__version__,source_page='https://snap.stanford.edu/data/email-Eu-core.html',source_url=url,sha256=hashlib.sha256(raw).hexdigest(),records=records,self_loop_records_removed=loops,nodes=len(g),simple_undirected_edges=g.number_of_edges(),transformation='Directed data projected to simple undirected graph; reciprocal records collapsed, loops removed, endpoints retained. This changes the estimand; no directed-network claim.',truth=target.tolist(),interpretation='Each stored graph is treated as a fixed artificial population. Sampling probabilities are assigned by the experimenter, not inferred from the real collection. No users, data-acquisition validation or interval coverage claim.',results=[])
for regime in ['uniform_0.2','heterogeneous_0.1_to_0.3']:
 p={v:(.2 if regime=='uniform_0.2' else .1+.2*(v%11)/10) for v in g}
 records=[]
 for seed in range(300):
  rng=np.random.default_rng(91000+seed);s=g.subgraph([v for v in g if rng.random()<p[v]])
  start=time.perf_counter();r=calc(s,p);duration=time.perf_counter()-start
  records.append(dict(seed=91000+seed,nodes=len(s),counts=r['observed_counts'],estimate=r['estimated_totals'],covariance=r['estimated_covariance'],seconds=duration,indefinite=any('not_positive_semidefinite' in w for w in r['warnings'])))
  if seed%50==49:print(regime,seed+1,'/ 300',flush=True)
 a=np.array([r['estimate'] for r in records]);c=np.array([r['covariance'] for r in records]);mu=a.mean(0);emp=np.cov(a,rowvar=False);se=a.std(0,ddof=1)/math.sqrt(len(a))
 # Bootstrap only quantifies Monte Carlo uncertainty of a variance ratio.
 # These are not confidence intervals for population motif totals.
 rng=np.random.default_rng(91919);ratios=[]
 for b in range(1000):
  ix=rng.integers(0,len(a),len(a));ratios.append(np.diag(c[ix].mean(0))/a[ix].var(0,ddof=1))
 result=dict(regime=regime,repetitions=300,mean_estimate=mu.tolist(),relative_bias=((mu-target)/target).tolist(),bias_in_mcse=((mu-target)/se).tolist(),empirical_covariance=emp.tolist(),mean_estimated_covariance=c.mean(0).tolist(),mean_variance_to_empirical_ratio=(np.diag(c.mean(0))/np.diag(emp)).tolist(),monte_carlo_variance_ratio_bootstrap_95=np.quantile(ratios,[.025,.975],axis=0).tolist(),indefinite_count=sum(r['indefinite'] for r in records),median_seconds=float(np.median([r['seconds'] for r in records])),raw=records)
 out['results'].append(result);(B/'public_study_results.json').write_text(json.dumps(out,indent=2,allow_nan=False))
 print(json.dumps({k:v for k,v in result.items() if k!='raw'}),flush=True)
