import json
from pathlib import Path
import numpy as np
B=Path(__file__).resolve().parent
actual=json.loads((B/'native_results.json').read_text())
expected={r['name']:r['result'] for r in json.loads((B/'native_expected.json').read_text())}
records=[]
for r in actual['results']:
 e=expected[r['name']];a=np.array(r['covariance']);b=np.array(e['estimated_covariance'])
 np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-8)
 np.testing.assert_allclose(r['totals'],e['estimated_totals'],rtol=1e-8,atol=1e-8)
 records.append(dict(name=r['name'],method=r['method'],maximum_scaled_covariance_error=float(np.max(abs(a-b)/np.maximum(1,abs(b)))),passed=True))
bench=[]
for name in expected:
 if name.startswith('correctness'):continue
 for method in ['sampling','survey']:
  rows=[r for r in actual['results'] if r['name']==name and r['method']==method]
  if rows:bench.append(dict(name=name,method=method,units=rows[0]['units'],repeats=len(rows),median_seconds=float(np.median([r['total_seconds'] for r in rows])),median_enumeration_seconds=float(np.median([r['enumeration_seconds'] for r in rows])),median_setup_seconds=float(np.median([r['setup_seconds'] for r in rows])),median_inference_seconds=float(np.median([r['inference_seconds'] for r in rows]))))
out=dict(R=actual['R'],versions=actual['versions'],comparisons=len(records),passed=len(records),distinct_graphs=len(expected),maximum_scaled_covariance_error=max(r['maximum_scaled_covariance_error'] for r in records),checks=records,benchmarks=bench,scope='Native installed sampling::varHT and survey::svytotal with sparse ppscov, wrapped in our original graph-to-motif adapter. Adapter and native stages separately timed. Process imports and JSON I/O excluded. R timings have coarse millisecond/centisecond resolution on this machine; zeros do not mean zero work.',harness_issue='First run stopped because the survey design API rejects a single observed motif unit. Correctness fixtures then required at least two edges; this is not a MotifSurvey failure. The package itself supports a single edge.')
(B/'native_summary.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ('checks',)}))
