"""Fresh native R process per fixture/method; actual Windows peak working set."""
import ctypes,json,subprocess,time,os,shutil
import psutil
from pathlib import Path
B=Path(__file__).resolve().parent
runtime=Path(os.environ.get('MOTIFSURVEY_RSCRIPT',str(B/'R-4.6.1/bin/Rscript.exe')))
if not runtime.exists():runtime=Path(shutil.which('Rscript') or '')
if not runtime.is_file():raise RuntimeError('Set MOTIFSURVEY_RSCRIPT to an installed Rscript executable')
library=os.environ.get('MOTIFSURVEY_RLIB',str(B/'r-library'))
class PMC(ctypes.Structure):
 _fields_=[('cb',ctypes.c_ulong),('PageFaultCount',ctypes.c_ulong)]+[(name,ctypes.c_size_t) for name in ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
ctypes.windll.psapi.GetProcessMemoryInfo.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_ulong]
fixtures=json.loads((B/'native_fixtures.json').read_text());records=[]
for fixture in fixtures:
 if fixture['name'].startswith('correctness'):continue
 for method in fixture['methods']:
  d=dict(fixture,methods=[method]);stem='native_isolated_'+fixture['name']+'_'+method
  source=B/(stem+'_input.json');dest=B/(stem+'_output.json');source.write_text(json.dumps([d]))
  command=[str(runtime),'--vanilla',str(B/'native_compare.R'),library,str(source),str(dest)]
  start=time.perf_counter();peak=0;process_peaks={}
  with (B/(stem+'.log')).open('w') as log:
   process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
   while process.poll() is None:
    try:
     parent=psutil.Process(process.pid)
     for proc in [parent]+parent.children(recursive=True):
      try:
       value=proc.memory_info().peak_wset;key=str(proc.pid)+':'+proc.name()
       process_peaks[key]=max(value,process_peaks.get(key,0));peak=max(peak,value)
      except psutil.Error:pass
    except psutil.Error:pass
    if time.perf_counter()-start>120:process.kill();raise TimeoutError(stem)
    time.sleep(.025)
  if process.returncode:raise RuntimeError((B/(stem+'.log')).read_text())
  result=json.loads(dest.read_text());r=dict(name=fixture['name'],method=method,largest_process_peak_working_set=peak,process_peaks=process_peaks,process_wall_seconds=time.perf_counter()-start,results=result['results'])
  assert peak>20_000_000, 'Only launcher observed: not a valid R memory measurement'
  records.append(r);(B/'native_memory_results.json').write_text(json.dumps(dict(measurement='psutil 7.2.2 Windows peak_wset of launcher and descendants, polled every 25ms. Report largest process peak, not a sum of asynchronous peaks. Includes runtime libraries and input. Timings within R exclude loading and JSON. Peak polling can miss final process activity.',corrected_harness='Initial attempt measured only the small Rscript launcher. Those numbers were rejected and all fixtures rerun with descendant process monitoring.',records=records),indent=2))
  print(json.dumps({k:v for k,v in r.items() if k!='results'}),flush=True)
