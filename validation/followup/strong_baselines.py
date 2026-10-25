"""New local SciPy adapters, not existing motif software or independent novelty."""
from itertools import combinations
from math import prod
import numpy as np
from scipy.sparse import coo_matrix,diags

def sparse_overlap(rows,p):
    """Generic sparse weighted sampling-covariance, as accepted by R survey.

    Cost depends on nonzero overlapping motif pairs; potentially quadratic.
    """
    m=len(rows);n=max(p,default=-1)+1
    if not m:return np.zeros(3),np.zeros((3,3))
    ii=[];jj=[];weights=[];types=[]
    for i,(s,t) in enumerate(rows):
        ii.extend([i]*len(s));jj.extend(s)
        weights.append(1/prod(p[v] for v in s));types.append(t)
    a=coo_matrix((np.ones(len(ii)),(ii,jj)),shape=(m,n)).tocsr()
    kernel=(a@diags(np.log([p[v] for v in range(n)]))@a.T).tocsr()
    kernel.data=-np.expm1(kernel.data)
    y=coo_matrix((weights,(range(m),types)),shape=(m,3)).tocsr()
    return np.asarray(y.sum(axis=0)).ravel(),(y.T@kernel@y).toarray()

def vectorized_subsets(rows,p):
    """Same identity as MotifSurvey; sparse duplicate-sum and BLAS reduction.

    This checks implementation overhead, not whether the method is already
    supplied by another released package. Enumeration still done by caller.
    """
    keys={};ii=[];jj=[];vv=[];coeff=[];total=np.zeros(3)
    for s,t in rows:
        weight=1/prod(p[v] for v in s);total[t]+=weight;s=tuple(sorted(s))
        for k in range(1,len(s)+1):
            for subset in combinations(s,k):
                ix=keys.get(subset)
                if ix is None:
                    ix=len(keys);keys[subset]=ix
                    coeff.append((-1)**(k+1)*prod(1-p[v] for v in subset))
                ii.append(ix);jj.append(t);vv.append(weight)
    if not keys:return total,np.zeros((3,3))
    g=coo_matrix((vv,(ii,jj)),shape=(len(keys),3)).toarray()
    return total,g.T@(g*np.array(coeff)[:,None])
