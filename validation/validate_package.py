import gzip,hashlib,json,math,platform,sys,time,urllib.request
from pathlib import Path
from itertools import combinations,product
import numpy as np
import networkx as nx
from motifsurvey import estimate
B=Path(__file__).resolve().parent;ROOT=B.parents[1]
def calc(g,p):
 return estimate([str(v) for v in g],[(str(u),str(v)) for u,v in g.edges],{str(v):p[v] for v in g},design='independent_node',complete_induced=True)
def truth(g):return np.array([g.number_of_edges(),sum(d*(d-1)//2 for _,d in g.degree()),sum(nx.triangles(g).values())//3],float)
def dump(name,value):(B/name).write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')
start=time.perf_counter();rng=np.random.default_rng(20261025);oracle=[]
for case in range(12):
 n=5+case%3;g=nx.Graph();g.add_nodes_from(range(n));g.add_edges_from(e for e in combinations(range(n),2) if rng.random()<.3+.1*(case%4))
 p={v:float(rng.uniform(.2,.9)) for v in g};mu=np.zeros(3);second=np.zeros((3,3));ecov=np.zeros((3,3));target=truth(g)
 for bits in product([0,1],repeat=n):
  probability=math.prod(p[v] if bits[v] else 1-p[v] for v in g);sample=g.subgraph([v for v in g if bits[v]])
  r=calc(sample,p);a=np.array(r['estimated_totals']);c=np.array(r['estimated_covariance']);mu+=probability*a;second+=probability*np.outer(a,a);ecov+=probability*c
 cov=second-np.outer(target,target)
 np.testing.assert_allclose(mu,target,rtol=1e-8,atol=1e-8);np.testing.assert_allclose(ecov,cov,rtol=1e-8,atol=1e-8)
 oracle.append(dict(case=case,samples=2**n,max_mean_error=float(abs(mu-target).max()),max_covariance_error=float(abs(ecov-cov).max())))
dump('package_oracle.json',oracle);print('Package oracle passed',flush=True)
# Exact counterexample to automatic normal confidence intervals; no such CI is in API.
g=nx.star_graph(5);p={v:.1 if v==0 else .5 for v in g};target=truth(g);coverage=np.zeros(3);zero=np.zeros(3)
for bits in product([0,1],repeat=6):
 probability=math.prod(p[v] if bits[v] else 1-p[v] for v in g);r=calc(g.subgraph([v for v in g if bits[v]]),p);a=np.array(r['estimated_totals']);d=np.diag(r['estimated_covariance']);sd=np.sqrt(np.maximum(d,0))
 coverage+=probability*(abs(a-target)<=1.95996398454*sd);zero+=probability*(a==0)
failure=dict(graph='5-leaf star',sampling='center p=.1, leaves p=.5',target=target.tolist(),naive_95pct_normal_interval_coverage=coverage.tolist(),probability_zero_observed_count=zero.tolist(),interpretation='Triangles truly absent; edge and wedge failure demonstrates that unbiased variance is not an interval guarantee. API returns no intervals.')
dump('normal_interval_failure.json',failure);print(json.dumps(failure),flush=True)
url='https://snap.stanford.edu/data/ca-GrQc.txt.gz';path=B/'ca-GrQc.txt.gz'
if not path.exists():
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'research-audit'}),timeout=45) as response:path.write_bytes(response.read())
raw=path.read_bytes();lines=gzip.decompress(raw).decode().splitlines();g=nx.Graph();records=0;loops=0
for line in lines:
 if not line or line.startswith('#'):continue
 u,v=map(int,line.split());records+=1;g.add_nodes_from((u,v))
 if u==v:loops+=1;continue
 g.add_edge(u,v)
target=truth(g);results=[];example_saved=False
for regime in ['uniform_0.3','heterogeneous_0.15_to_0.6']:
 p={v:(.3 if regime=='uniform_0.3' else .15+.45*(v%11)/10) for v in g}
 estimates=[];variances=[];durations=[];psd=[];counts=[]
 for seed in range(30):
  rng=np.random.default_rng(6000+seed);s=g.subgraph([v for v in g if rng.random()<p[v]])
  t=time.perf_counter();r=calc(s,p);durations.append(time.perf_counter()-t);estimates.append(r['estimated_totals']);variances.append(np.diag(r['estimated_covariance']));psd.append(r['diagnostics']['minimum_covariance_eigenvalue']>=-1e-8);counts.append(r['observed_counts'])
  if not example_saved:
   dump('public_example_result.json',r);example_saved=True
 a=np.array(estimates);v=np.array(variances);mcse=a.std(axis=0,ddof=1)/math.sqrt(len(a));bias=a.mean(axis=0)-target
 results.append(dict(regime=regime,repeats=30,mean_estimate=a.mean(axis=0).tolist(),relative_bias=(bias/target).tolist(),mean_monte_carlo_standard_error=mcse.tolist(),bias_in_mcse=(bias/mcse).tolist(),empirical_sampling_variance=a.var(axis=0,ddof=1).tolist(),mean_estimated_variance=v.mean(axis=0).tolist(),median_seconds=float(np.median(durations)),max_seconds=max(durations),indefinite_covariance_samples=30-sum(psd),raw_estimates=estimates,raw_estimated_variances=[x.tolist() for x in variances],observed_counts=counts))
 print(json.dumps({k:v for k,v in results[-1].items() if not k.startswith('raw') and k!='observed_counts'}),flush=True)
 dump('public_network_results.json',dict(source_url=url,source_page='https://snap.stanford.edu/data/ca-GrQc.html',sha256=hashlib.sha256(raw).hexdigest(),records=records,self_loop_records_removed=loops,nodes=g.number_of_nodes(),simple_undirected_edges=g.number_of_edges(),deduplication='Reciprocal and repeated edge records collapse to a single undirected edge; loops removed explicitly',population_truth=target.tolist(),design='Synthetic independent node thinning of a public graph; not evidence that original collection was random sampling',results=results))
dump('package_environment.json',dict(python=sys.version,numpy=np.__version__,networkx=nx.__version__,platform=platform.platform(),total_wall_seconds=time.perf_counter()-start))
