"""Independent-node design covariance feasibility test. All code original."""
from itertools import combinations,product
from collections import defaultdict
from math import prod,fsum
from pathlib import Path
import json,time,tracemalloc,platform
import numpy as np
B=Path(__file__).resolve().parent

def instances(nodes,edges):
    adj={v:set() for v in nodes}
    for u,v in edges:adj[u].add(v);adj[v].add(u)
    rows=[]
    for u,v in edges:rows.append((frozenset((u,v)),0))
    for u in nodes:
        for v,w in combinations(sorted(adj[u]),2):rows.append((frozenset((u,v,w)),1))
    for u in nodes:
        for v in adj[u]:
            if v>u:
                for w in adj[u]&adj[v]:
                    if w>v:rows.append((frozenset((u,v,w)),2))
    return rows

def aggregate(rows,p):
    sums={};total=np.zeros(3)
    for vertices,t in rows:
        weight=1/prod(p[v] for v in vertices);total[t]+=weight
        ordered=sorted(vertices)
        for k in range(1,len(ordered)+1):
            for subset in combinations(ordered,k):
                if subset not in sums:sums[subset]=np.zeros(3)
                sums[subset][t]+=weight
    cov=np.zeros((3,3))
    for subset,g in sums.items():
        coeff=(-1)**(len(subset)+1)*prod(1-p[v] for v in subset)
        cov+=coeff*np.outer(g,g)
    return total,cov

def direct(rows,p):
    total=np.zeros(3);cov=np.zeros((3,3))
    weights=[1/prod(p[v] for v in s) for s,t in rows]
    for i,(s,t) in enumerate(rows):
        total[t]+=weights[i]
        for j,(a,b) in enumerate(rows):
            cov[t,b]+=weights[i]*weights[j]*(1-prod(p[v] for v in s&a))
    return total,cov

def indexed(rows,p):
    bynode=defaultdict(list);total=np.zeros(3);cov=np.zeros((3,3))
    weights=[1/prod(p[v] for v in s) for s,t in rows]
    for i,(s,t) in enumerate(rows):
        total[t]+=weights[i]
        for v in s:bynode[v].append(i)
    for i,(s,t) in enumerate(rows):
        js=set().union(*(bynode[v] for v in s))
        for j in js:
            a,b=rows[j];cov[t,b]+=weights[i]*weights[j]*(1-prod(p[v] for v in s&a))
    return total,cov

def vectorized(rows,p):
    # Dense matrix baseline, avoids Python pair-loop overhead. Matrix-building included.
    n=max(p,default=-1)+1;m=len(rows)
    incidence=np.zeros((m,n));y=np.zeros((m,3))
    for i,(s,t) in enumerate(rows):incidence[i,list(s)]=1;y[i,t]=1/prod(p[v] for v in s)
    kernel=-np.expm1((incidence*np.log([p[v] for v in range(n)]))@incidence.T)
    return y.sum(axis=0),y.T@kernel@y

def check():
    rng=np.random.default_rng(20261025);exact=[];agreements=0
    for case in range(12):
        n=5+case%3;nodes=list(range(n));edges=[(u,v) for u,v in combinations(nodes,2) if rng.random()<.3+.1*(case%4)]
        p={v:float(rng.uniform(.2,.9)) for v in nodes};allrows=instances(nodes,edges)
        truth=np.bincount([t for s,t in allrows],minlength=3).astype(float)
        mean=np.zeros(3);second=np.zeros((3,3));mean_cov=np.zeros((3,3));samples=0
        for bits in product((0,1),repeat=n):
            keep={v for v in nodes if bits[v]};prob=prod(p[v] if bits[v] else 1-p[v] for v in nodes)
            rows=instances(sorted(keep),[(u,v) for u,v in edges if u in keep and v in keep])
            a,c=aggregate(rows,p);b,d=direct(rows,p)
            assert np.allclose(a,b,rtol=1e-8,atol=1e-8) and np.allclose(c,d,rtol=1e-8,atol=1e-8)
            mean+=prob*a;second+=prob*np.outer(a,a);mean_cov+=prob*c;samples+=1;agreements+=1
        true_cov=second-np.outer(truth,truth)
        assert np.allclose(mean,truth,rtol=1e-8,atol=1e-8) and np.allclose(mean_cov,true_cov,rtol=1e-8,atol=1e-8)
        exact.append(dict(case=case,nodes=n,edges=len(edges),samples=samples,mean_error=float(np.max(abs(mean-truth))),covariance_error=float(np.max(abs(mean_cov-true_cov)))))
    for case in range(30):
        nodes=list(range(10));edges=[e for e in combinations(nodes,2) if rng.random()<.4];p={v:float(rng.uniform(.1,1)) for v in nodes};rows=instances(nodes,edges)
        a,c=aggregate(rows,p)
        for fn in (direct,indexed,vectorized):
            b,d=fn(rows,p);assert np.allclose(a,b,rtol=1e-8,atol=1e-8) and np.allclose(c,d,rtol=1e-8,atol=1e-8);agreements+=1
    return dict(exhaustive=exact,baseline_agreements=agreements)

def measure(fn,rows,p,repeats=3):
    timings=[]
    for _ in range(repeats):
        t=time.perf_counter();result=fn(rows,p);timings.append(time.perf_counter()-t)
    tracemalloc.start();fn(rows,p);_,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
    return dict(seconds=timings,median_seconds=float(np.median(timings)),peak_traced_bytes=peak,covariance=result[1].tolist(),total=result[0].tolist())

if __name__=='__main__':
    result=check();(B/'preflight_correctness.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
    rng=np.random.default_rng(440);bench=[]
    fixtures=[('ring',2000,[(i,(i+1)%2000) for i in range(2000)]),('hub',1000,[(0,i) for i in range(1,1000)]),('random',700,[(i,j) for i,j in combinations(range(700),2) if rng.random()<.004])]
    for name,n,edges in fixtures:
        p={i:float(.3+.6*(i%7)/6) for i in range(n)};rows=instances(list(range(n)),edges)
        # Hub has ~500k wedges: benchmark covariance on edge units, then all motifs for fast path.
        # Full hub pair baseline explicitly omitted at this size, not reported as an executed speedup.
        if name=='hub':compared=[r for r in rows if r[1]==0]
        else:compared=rows
        entry=dict(name=name,nodes=n,edges=len(edges),all_motif_instances=len(rows),compared_instances=len(compared),units='edges only' if name=='hub' else 'edges, wedges, triangles')
        for fn in (aggregate,indexed,vectorized):
            entry[fn.__name__]=measure(fn,compared,p)
        for fn in ('indexed','vectorized'):
            assert np.allclose(entry['aggregate']['covariance'],entry[fn]['covariance'],rtol=1e-8,atol=1e-8)
        bench.append(entry);(B/'preflight_benchmark.json').write_text(json.dumps(bench,indent=2));print(json.dumps(entry),flush=True)
    (B/'preflight_environment.json').write_text(json.dumps(dict(python=platform.python_version(),platform=platform.platform(),numpy=np.__version__,baseline='Local Python implementations, not native R sampling execution',memory='tracemalloc, not OS peak RSS'),indent=2))
