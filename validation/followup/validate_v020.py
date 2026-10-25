import gzip,hashlib,inspect,json,math,sys,urllib.request,subprocess,re,time
from itertools import product
from pathlib import Path
import networkx as nx
import numpy as np
import motifsurvey
import motifsurvey.core
from stress import calc
B=Path(__file__).resolve().parent
results=dict(version=motifsurvey.__version__,native_fixture_checks=[],exact_samples=0,public_replay_checks=0)
native=json.loads((B/'native_results.json').read_text())
fixtures=json.loads((B/'native_fixtures.json').read_text())
for f in fixtures:
 g=nx.Graph();g.add_nodes_from(range(f['n']));g.add_edges_from(f['edges']);p=dict(enumerate(f['p']));r=calc(g,p)
 for method in f['methods']:
  n=next(x for x in native['results'] if x['name']==f['name'] and x['method']==method)
  np.testing.assert_allclose(r['estimated_totals'],n['totals'],rtol=1e-8,atol=1e-8)
  np.testing.assert_allclose(r['estimated_covariance'],n['covariance'],rtol=1e-8,atol=1e-8)
  results['native_fixture_checks'].append(dict(name=f['name'],method=method,passed=True))
# New reduction code: repeat the design expectation oracle, without rerunning
# the original public experiment or native binaries whose inputs are unchanged.
for case in range(12):
 n=5+case%3;g=nx.fast_gnp_random_graph(n,.2+.15*(case%4),seed=5500+case)
 p={v:.2+.65*((v+case)%7)/6 for v in g}
 truth=np.array([g.number_of_edges(),sum(d*(d-1)//2 for _,d in g.degree()),sum(nx.triangles(g).values())//3],float)
 mu=np.zeros(3);second=np.zeros((3,3));ecov=np.zeros((3,3))
 for bits in product([0,1],repeat=n):
  probability=math.prod(p[v] if bits[v] else 1-p[v] for v in g);r=calc(g.subgraph([v for v in g if bits[v]]),p)
  a=np.array(r['estimated_totals']);c=np.array(r['estimated_covariance']);mu+=probability*a;second+=probability*np.outer(a,a);ecov+=probability*c;results['exact_samples']+=1
 np.testing.assert_allclose(mu,truth,rtol=1e-8,atol=1e-8)
 np.testing.assert_allclose(ecov,second-np.outer(truth,truth),rtol=1e-8,atol=1e-8)
study=json.loads((B/'public_study_results.json').read_text());g=nx.Graph()
if not (B/'email-Eu-core.txt.gz').exists():
 with urllib.request.urlopen(study['source_url'],timeout=45) as response:(B/'email-Eu-core.txt.gz').write_bytes(response.read())
assert hashlib.sha256((B/'email-Eu-core.txt.gz').read_bytes()).hexdigest()==study['sha256']
for line in gzip.decompress((B/'email-Eu-core.txt.gz').read_bytes()).decode().splitlines():
 if not line or line.startswith('#'):continue
 u,v=map(int,line.split());g.add_nodes_from((u,v))
 if u!=v:g.add_edge(u,v)
for regime in study['results']:
 p={v:(.2 if regime['regime']=='uniform_0.2' else .1+.2*(v%11)/10) for v in g}
 for record in regime['raw'][::10]:
  rng=np.random.default_rng(record['seed']);s=g.subgraph([v for v in g if rng.random()<p[v]]);r=calc(s,p)
  np.testing.assert_allclose(r['estimated_totals'],record['estimate'],rtol=1e-8,atol=1e-8)
  np.testing.assert_allclose(r['estimated_covariance'],record['covariance'],rtol=1e-8,atol=1e-8);results['public_replay_checks']+=1
source=Path(inspect.getfile(motifsurvey.core));delivered=B.parents[1]/'prototypes/motifsurvey/src/motifsurvey/core.py'
if not delivered.exists():delivered=B.parents[1]/'src/motifsurvey/core.py'
results['installed_core_sha256']=hashlib.sha256(source.read_bytes()).hexdigest();results['delivered_core_sha256']=hashlib.sha256(delivered.read_bytes()).hexdigest()
assert results['installed_core_sha256']==results['delivered_core_sha256']
test_start=time.perf_counter()
test=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-v'],cwd=delivered.parents[2],capture_output=True,text=True)
if test.returncode:raise RuntimeError(test.stderr)
count=int(re.search(r'Ran (\d+) tests',test.stderr).group(1))
results['unit_tests']=dict(tests=count,passed=count,wall_seconds=time.perf_counter()-test_start,command='python -m unittest discover -s tests -v',output=test.stderr)
results['all_passed']=True
(B/'v020_validation.json').write_text(json.dumps(results,indent=2))
print(json.dumps({k:v for k,v in results.items() if k!='native_fixture_checks'}),flush=True)
