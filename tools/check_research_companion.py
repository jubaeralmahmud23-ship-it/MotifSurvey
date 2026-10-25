"""Record the four companion tests separately from frozen release checks."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys

import motifsurvey
import motifsurvey.core

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
args.output.mkdir(parents=True, exist_ok=True)
process = subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "validation/research_usefulness", "-p", "test_*.py", "-v"], cwd=root, capture_output=True, text=True, encoding="utf-8")
log = process.stdout + process.stderr
print(log)
count = re.search(r"Ran (\d+) tests", log)
report = {"passed": process.returncode == 0 and count is not None and int(count.group(1)) == 4,
          "tests": int(count.group(1)) if count else None, "python": platform.python_version(), "system": platform.system(),
          "software_version": motifsurvey.__version__, "github_sha": os.environ.get("GITHUB_SHA"),
          "core_sha256": hashlib.sha256(Path(motifsurvey.core.__file__).read_bytes()).hexdigest(),
          "scope": "Additional companion tests; cached native outputs are checked, R is not executed here"}
(args.output / "companion.log").write_text(log, encoding="utf-8")
(args.output / "companion.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
sys.exit(0 if report["passed"] else 1)
