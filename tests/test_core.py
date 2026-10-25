import copy,itertools,json,subprocess,sys,tempfile,unittest
from pathlib import Path
import numpy as np
from motifsurvey import estimate, ResourceLimitError

def payload():return dict(nodes=['a','b','c'],edges=[['a','b'],['b','c'],['a','c']],probabilities=dict(a=.5,b=.5,c=.5),design='independent_node',complete_induced=True)

class TestMotifSurvey(unittest.TestCase):
    def test_single_edge_closed_form(self):
        r=estimate(['a','b'],[('a','b')],dict(a=.5,b=.25),design='independent_node',complete_induced=True)
        self.assertEqual(r['estimated_totals'],[8,0,0]);self.assertEqual(r['estimated_covariance'][0][0],56)
    def test_triangle_multiplicities(self):
        r=estimate(**payload());self.assertEqual(r['observed_counts'],[3,3,1]);self.assertEqual(r['estimated_totals'],[12,24,8])
        self.assertAlmostEqual(r['estimated_covariance'][1][1],504);self.assertAlmostEqual(r['estimated_covariance'][2][2],56);self.assertAlmostEqual(r['estimated_covariance'][1][2],168)
    def test_census(self):
        x=payload();x['probabilities']=dict(a=1,b=1,c=1);r=estimate(**x)
        self.assertEqual(r['estimated_totals'],[3,3,1]);np.testing.assert_array_equal(r['estimated_covariance'],np.zeros((3,3)))
    def test_empty_sample(self):
        r=estimate([],[],{},design='independent_node',complete_induced=True);self.assertEqual(r['estimated_totals'],[0,0,0]);self.assertTrue(r['warnings'])
    def test_isolates(self):
        r=estimate(['a'],[],{'a':.5},design='independent_node',complete_induced=True);self.assertEqual(r['diagnostics']['observed_nodes'],1)
    def test_permutation_invariance(self):
        x=payload();a=estimate(**x);x['nodes'].reverse();x['edges']=[e[::-1] for e in reversed(x['edges'])];b=estimate(**x);self.assertEqual(a,b)
    def test_wrong_designs(self):
        for design in ['fixed_size','ego','snowball','edge','unknown',None]:
            x=payload();x['design']=design
            with self.subTest(design=design),self.assertRaises(ValueError):estimate(**x)
    def test_incomplete_edges(self):
        x=payload();x['complete_induced']=False
        with self.assertRaises(ValueError):estimate(**x)
    def test_bad_probabilities(self):
        for p in [0,-1,1.1,float('inf'),float('nan'),True,'0.5',None]:
            x=payload();x['probabilities']['a']=p
            with self.subTest(p=p),self.assertRaises(ValueError):estimate(**x)
    def test_probability_keys(self):
        for p in [dict(a=.5,b=.5),dict(a=.5,b=.5,c=.5,d=.5)]:
            x=payload();x['probabilities']=p
            with self.assertRaises(ValueError):estimate(**x)
    def test_rejects_nonreal_probability_scalars(self):
        for value in [np.bool_(True),np.bool_(False),complex(.5,0),np.complex128(.5+.2j),np.complex64(.5),b'0.5',bytearray(b'0.5'),np.array(.5)]:
            x=payload();x['probabilities']['a']=value
            with self.subTest(value_type=type(value).__name__),self.assertRaises(ValueError):estimate(**x)
    def test_huge_numeric_probability_is_clean_error(self):
        x=payload();x['probabilities']['a']=10**5000
        with self.assertRaises(ValueError):estimate(**x)
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'input.json'
            # 401 digits are below Python's JSON digit limit but overflow float.
            x['probabilities']['a']=10**400;source.write_text(json.dumps(x))
            process=subprocess.run([sys.executable,'-m','motifsurvey',str(source)],capture_output=True,text=True)
            self.assertEqual(process.returncode,2)
            self.assertNotIn('Traceback',process.stderr)
    def test_supported_real_scalar_probabilities(self):
        from fractions import Fraction
        from decimal import Decimal
        expected=estimate(**payload())['estimated_totals']
        for value in [np.float32(.5),np.float64(.5),Decimal('0.5'),Fraction(1,2)]:
            x=payload();x['probabilities']['a']=value
            self.assertEqual(estimate(**x)['estimated_totals'],expected)
    def test_duplicate_nodes(self):
        x=payload();x['nodes'].append('a')
        with self.assertRaises(ValueError):estimate(**x)
    def test_invalid_ids(self):
        for nodes in [[''],[1],[None]]:
            with self.assertRaises(ValueError):estimate(nodes,[],{},design='independent_node',complete_induced=True)
    def test_loops(self):
        x=payload();x['edges'].append(['a','a'])
        with self.assertRaises(ValueError):estimate(**x)
    def test_duplicate_edges(self):
        x=payload();x['edges'].append(['b','a'])
        with self.assertRaises(ValueError):estimate(**x)
    def test_unknown_endpoint(self):
        x=payload();x['edges'].append(['a','d'])
        with self.assertRaises(ValueError):estimate(**x)
    def test_bad_edge_shapes(self):
        for e in [['a'],['a','b','c'],'ab']:
            x=payload();x['edges']=[e]
            with self.assertRaises(ValueError):estimate(**x)
    def test_resource_budgets(self):
        with self.assertRaises(ResourceLimitError):estimate(**payload(),max_instances=6)
        with self.assertRaises(ResourceLimitError):estimate(**payload(),max_subsets=1)
        with self.assertRaises(ValueError):estimate(**payload(),max_instances=True)
    def test_underflow(self):
        x=payload();x['probabilities']=dict(a=1e-200,b=1e-200,c=1e-200)
        with self.assertRaises(ValueError):estimate(**x)
    def test_indefinite_covariance_retained(self):
        n=20;nodes=list(map(str,range(n)));r=estimate(nodes,[(str(i),str((i+1)%n)) for i in range(n)],dict.fromkeys(nodes,.5),design='independent_node',complete_induced=True)
        self.assertLess(r['diagnostics']['minimum_covariance_eigenvalue'],0)
        self.assertTrue(any('not_positive_semidefinite' in w for w in r['warnings']))
    def test_strict_json(self):json.dumps(estimate(**payload()),allow_nan=False)
    def test_multiple_blocks_against_pair_formula(self):
        # This graph crosses several reduction blocks. The oracle explicitly
        # sums ordered motif pairs rather than repeating subset aggregation.
        n=300;nodes=list(map(str,range(n)))
        edges=[(str(i),str((i+1)%n)) for i in range(n)]
        p={v:.2+.7*(i%7)/6 for i,v in enumerate(nodes)}
        r=estimate(nodes,edges,p,design='independent_node',complete_induced=True)
        self.assertGreater(r['diagnostics']['subset_keys'],1024)
        rows=[(frozenset(e),0) for e in edges]+[(frozenset((str((i-1)%n),str(i),str((i+1)%n))),1) for i in range(n)]
        from math import prod,fsum
        weights=[1/prod(p[v] for v in s) for s,t in rows]
        oracle=np.zeros((3,3))
        for t in range(2):
            for u in range(2):
                oracle[t,u]=fsum(weights[i]*weights[j]*(1-prod(p[v] for v in s&a)) for i,(s,ty) in enumerate(rows) if ty==t for j,(a,uy) in enumerate(rows) if uy==u)
        np.testing.assert_allclose(r['estimated_covariance'],oracle,rtol=1e-12,atol=1e-9)
        np.testing.assert_array_equal(r['estimated_covariance'],np.asarray(r['estimated_covariance']).T)
    def test_cli_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            src=Path(d)/'input.json';out=Path(d)/'output.json';src.write_text(json.dumps(payload()))
            cmd=[sys.executable,'-m','motifsurvey',str(src),'--output',str(out)]
            a=subprocess.run(cmd,capture_output=True,text=True);self.assertEqual(a.returncode,0,a.stderr)
            before=out.read_bytes();b=subprocess.run(cmd,capture_output=True,text=True);self.assertEqual(b.returncode,2);self.assertEqual(out.read_bytes(),before)
    def test_networkx_native_counts_all_4_node_graphs(self):
        import networkx as nx
        vertices=list(map(str,range(4)));possible=list(itertools.combinations(vertices,2))
        for bits in itertools.product([0,1],repeat=6):
            g=nx.Graph();g.add_nodes_from(vertices);g.add_edges_from(e for e,b in zip(possible,bits) if b)
            r=estimate(vertices,g.edges,dict.fromkeys(vertices,1),design='independent_node',complete_induced=True)
            expected=[g.number_of_edges(),sum(d*(d-1)//2 for _,d in g.degree()),sum(nx.triangles(g).values())//3]
            self.assertEqual(r['observed_counts'],expected)

if __name__=='__main__':unittest.main()
