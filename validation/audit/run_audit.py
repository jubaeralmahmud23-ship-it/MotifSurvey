"""Independent rational-arithmetic and reported-result audit; MIT License."""
import argparse
from fractions import Fraction as F
from itertools import combinations, product
from math import prod
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import motifsurvey
import motifsurvey.core

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'examples'))
from closure_profile import closure_profile


def motifs(vertices, edges):
    edges={frozenset(e) for e in edges}
    units=[(e,0) for e in edges]
    # Enumerate triples and their possible centers, independently of the core.
    for triple in combinations(vertices,3):
        for center in triple:
            others=[v for v in triple if v!=center]
            if all(frozenset((center,v)) in edges for v in others):
                units.append((frozenset(triple),1))
        if all(frozenset(e) in edges for e in combinations(triple,2)):
            units.append((frozenset(triple),2))
    return units


def rational_audit():
    vertices=list('abcd'); possible=list(combinations(vertices,2))
    regimes=[[F(1,2)]*4,[F(1,4),F(1,2),F(3,4),F(1)]]
    checked=0;maximum=0.
    for mask in range(64):
        edges=[e for i,e in enumerate(possible) if mask&(1<<i)]
        units=motifs(vertices,edges)
        truth=[sum(t==k for _,t in units) for k in range(3)]
        for probabilities in regimes:
            p=dict(zip(vertices,probabilities))
            weights=[1/prod(p[v] for v in s) for s,t in units]
            mean=[F(0)]*3;second=[[F(0) for _ in range(3)] for _ in range(3)]
            mean_cov=[[F(0) for _ in range(3)] for _ in range(3)]
            for bits in product((0,1),repeat=4):
                selected={v for v,b in zip(vertices,bits) if b}
                chance=prod(p[v] if v in selected else 1-p[v] for v in vertices)
                indices=[i for i,(s,t) in enumerate(units) if s<=selected]
                totals=[sum((weights[i] for i in indices if units[i][1]==k),F(0)) for k in range(3)]
                covariance=[[F(0) for _ in range(3)] for _ in range(3)]
                for i in indices:
                    s,t=units[i]
                    for j in indices:
                        a,u=units[j]
                        covariance[t][u]+=weights[i]*weights[j]*(1-prod(p[v] for v in s&a))
                result=motifsurvey.estimate(sorted(selected),[e for e in edges if set(e)<=selected],{v:float(p[v]) for v in selected},design='independent_node',complete_induced=True)
                a=np.array(result['estimated_covariance']);wanted=np.array([[float(x) for x in row] for row in covariance])
                maximum=max(maximum,float(np.max(abs(a-wanted)/np.maximum(1,abs(wanted)))))
                np.testing.assert_allclose(result['estimated_totals'],list(map(float,totals)),rtol=1e-12,atol=1e-10)
                np.testing.assert_allclose(a,wanted,rtol=1e-12,atol=1e-10)
                for t in range(3):
                    mean[t]+=chance*totals[t]
                    for u in range(3):
                        second[t][u]+=chance*totals[t]*totals[u]
                        mean_cov[t][u]+=chance*covariance[t][u]
                checked+=1
            assert mean==truth
            assert all(mean_cov[t][u]==second[t][u]-truth[t]*truth[u] for t in range(3) for u in range(3))
    return {'graphs':64,'probability_regimes':2,'exact_sampling_laws':128,'sample_calculations':checked,
            'expectation_identities_exact_in_fraction_arithmetic':True,'maximum_scaled_sample_covariance_error':maximum}


def reported_audit():
    folder=ROOT/'validation/design_case/results'
    raw=json.loads((folder/'realizations.json').read_text())
    summary=json.loads((folder/'summary.json').read_text())
    c=np.array([0,-.25,3.]);target=summary['metadata']['population']['contrast']
    maximum=0.;checks=0
    pilot={};mses={}
    for stage,designs in raw.items():
        for name,rows in designs.items():
            x=np.array([r['estimated_totals'] for r in rows]);v=np.array([r['estimated_covariance'] for r in rows])
            contrasts=x@c;variances=np.einsum('i,nij,j->n',c,v,c)
            np.testing.assert_allclose(contrasts,[r['contrast']['estimate'] for r in rows],rtol=1e-12,atol=1e-7)
            np.testing.assert_allclose(variances,[r['contrast']['variance'] for r in rows],rtol=1e-11,atol=1e-4)
            if stage=='pilot':
                pilot[name]=float(variances.mean())
            else:
                expected=summary['evaluation'][name]
                mses[name]=float(np.mean((contrasts-target)**2))
                np.testing.assert_allclose(mses[name],expected['empirical_mse'],rtol=1e-12)
                np.testing.assert_allclose(variances.mean(),expected['mean_estimated_variance'],rtol=1e-12)
                assert int((variances<0).sum())==expected['negative_estimated_variance']
                assert sum(any('not_positive_semidefinite' in w for w in r['warnings']) for r in rows)==expected['indefinite_covariance']
                exact=expected['exact_variance']
                mcse=variances.std(ddof=1)/np.sqrt(len(variances))
                np.testing.assert_allclose((variances.mean()-exact)/mcse,expected['estimated_variance_mean_error_in_mcse'],rtol=1e-12)
                maximum=max(maximum,abs(float((variances.mean()-exact)/mcse)))
            checks+=len(rows)
    assert min(pilot,key=pilot.get)==summary['pilot_selection']['full_covariance_choice']
    for name,value in mses.items():
        np.testing.assert_allclose(value/mses['uniform'],summary['evaluation'][name]['mse_ratio_to_uniform'],rtol=1e-12)
    return {'recorded_realizations_checked':checks,'mse_and_variance_summaries_reproduced':True,
            'maximum_mean_variance_discrepancy_in_mcse':maximum,'scope':'Algebraic recomputation from all saved totals/covariances, not new sampling runs'}


def examples():
    triangle=motifsurvey.estimate(list('abc'),[('a','b'),('b','c'),('a','c')],dict.fromkeys('abc',.5),design='independent_node',complete_induced=True)
    np.testing.assert_allclose(triangle['estimated_covariance'],[[84,216,72],[216,504,168],[72,168,56]],rtol=1e-12)
    profile=closure_profile(triangle,[.25])['profile'][0]
    assert profile['estimate']==18 and profile['estimated_variance']==283.5 and profile['diagonal_only_variance']==535.5
    single=motifsurvey.estimate(list('ab'),[('a','b')],{'a':.5,'b':.25},design='independent_node',complete_induced=True)
    assert single['estimated_totals'][0]==8 and single['estimated_covariance'][0][0]==56
    return {'single_edge_unequal':single,'observed_triangle_half_probabilities':triangle,'closure_contrast':profile,
            'triangle_population_distribution':{'probability_observed':.125,'ht_triangle_when_observed':8,'ht_triangle_otherwise':0,'expectation':1,'true_variance':7,'estimated_variance_when_observed':56,'expected_estimated_variance':7}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result={'version':motifsurvey.__version__,'core_sha256':hashlib.sha256(Path(motifsurvey.core.__file__).read_bytes()).hexdigest(),
            'rational_audit':rational_audit(),'reported_result_audit':reported_audit(),'teaching_examples':examples(),'passed':True}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in result.items() if k!='teaching_examples'},indent=2))
