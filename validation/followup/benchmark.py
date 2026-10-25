import argparse,json,os,platform,subprocess,sys,time,ctypes
from pathlib import Path
import numpy as np
import scipy
import networkx as nx
from motifsurvey import estimate
B=Path(__file__).resolve().parent
sys.path.insert(0,str(B.parent/'motif_search' if (B.parent/'motif_search').exists() else B.parent))
from preflight import instances
from strong_baselines import sparse_overlap,vectorized_subsets

def peak_memory():
 class PMC(ctypes.Structure):
  _fields_=[('cb',ctypes.c_ulong),('PageFaultCount',ctypes.c_ulong)]+[(name,ctypes.c_size_t) for name in ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
 x=PMC();x.cb=ctypes.sizeof(x)
 ctypes.windll.kernel32.GetCurrentProcess.restype=ctypes.c_void_p
 ctypes.windll.psapi.GetProcessMemoryInfo.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_ulong]
 handle=ctypes.windll.kernel32.GetCurrentProcess()
 if not ctypes.windll.psapi.GetProcessMemoryInfo(handle,ctypes.byref(x),x.cb):raise ctypes.WinError()
 return x.PeakWorkingSetSize

def fixture(name):
 if name=='ring_2000':return nx.cycle_graph(2000)
 if name=='hub_50':return nx.star_graph(49)
 if name=='hub_200':return nx.star_graph(199)
 if name=='hub_1000':return nx.star_graph(999)
 if name=='clique_25':return nx.complete_graph(25)
 if name=='heavy_tail_1500':return nx.barabasi_albert_graph(1500,3,seed=20261025)
 if name=='random_700':return nx.fast_gnp_random_graph(700,.004,seed=440)
 if name=='ring_10000':return nx.cycle_graph(10000)
 raise ValueError(name)

def worker(name,method):
 g=fixture(name);p={v:.3+.6*(v%7)/6 for v in g}
 before=peak_memory();times=[]
 def run():
  if method=='motifsurvey':
   r=estimate([str(v) for v in g],[(str(u),str(v)) for u,v in g.edges],{str(v):p[v] for v in g},design='independent_node',complete_induced=True,max_instances=2_000_000,max_subsets=4_000_000)
   return np.array(r['estimated_totals']),np.array(r['estimated_covariance'])
  rows=instances(list(g),list(g.edges));return {'sparse_overlap':sparse_overlap,'vectorized_subsets':vectorized_subsets}[method](rows,p)
 for _ in range(3):
  start=time.perf_counter();a,c=run();times.append(time.perf_counter()-start)
 return dict(fixture=name,method=method,nodes=len(g),edges=g.number_of_edges(),units=g.number_of_edges()+sum(d*(d-1)//2 for _,d in g.degree())+sum(nx.triangles(g).values())//3,seconds=times,median_seconds=float(np.median(times)),process_peak_working_set=peak_memory(),premeasurement_peak_working_set=before,totals=a.tolist(),covariance=c.tolist())

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--fixture');parser.add_argument('--method');args=parser.parse_args()
 if args.fixture:print(json.dumps(worker(args.fixture,args.method),allow_nan=False));sys.exit()
 out=dict(protocol='Fresh subprocess per fixture/method, three internal repetitions including motif enumeration. Library imports and process startup excluded from timing, included in peak process working set. All Python workers import identical libraries. Package includes validation/hash; local adapters do not. Limits raised explicitly to 2 million motif units/4 million subset keys.',versions=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,networkx=nx.__version__,platform=platform.platform()),records=[],skipped=[])
 for name in ['ring_2000','hub_50','random_700','clique_25','hub_200','heavy_tail_1500','ring_10000']:
  refs=[]
  for method in ['motifsurvey','sparse_overlap','vectorized_subsets']:
   # Potential huge M^2 overlap matrix: avoid uncontrolled multi-GiB allocation.
   if method=='sparse_overlap' and name in ('hub_200','heavy_tail_1500'):
    out['skipped'].append(dict(fixture=name,method=method,reason='Potential quadratic motif-overlap storage; not timed or reported as a measured speedup'));continue
   result=subprocess.run([sys.executable,__file__,'--fixture',name,'--method',method],capture_output=True,text=True,timeout=180)
   if result.returncode:raise RuntimeError(result.stderr)
   r=json.loads(result.stdout);out['records'].append(r);refs.append(r);print(json.dumps({k:v for k,v in r.items() if k not in ('totals','covariance')}),flush=True)
  for r in refs[1:]:
   np.testing.assert_allclose(r['totals'],refs[0]['totals'],rtol=1e-8,atol=1e-8)
   np.testing.assert_allclose(r['covariance'],refs[0]['covariance'],rtol=1e-8,atol=1e-8)
  (B/'strong_benchmark.json').write_text(json.dumps(out,indent=2,allow_nan=False))
