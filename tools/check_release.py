"""Portable checks of an installed distribution, with machine-readable evidence."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cache", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    # Keep small validation matrices from creating large BLAS thread pools.
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[variable] = "1"
    import motifsurvey
    import motifsurvey.core
    import networkx
    import numpy
    installed = Path(motifsurvey.core.__file__).resolve()
    if root / "src" in installed.parents:
        raise RuntimeError("Install the wheel or sdist; source/editable imports are not a release check")
    core_hash = hashlib.sha256(installed.read_bytes()).hexdigest()
    if core_hash != hashlib.sha256((root / "src/motifsurvey/core.py").read_bytes()).hexdigest():
        raise RuntimeError("Installed core differs from the delivered source")
    report = {
        "version": motifsurvey.__version__, "core_sha256": core_hash,
        "platform": platform.platform(), "system": platform.system(),
        "machine": platform.machine(), "python": platform.python_version(),
        "numpy": numpy.__version__, "networkx": networkx.__version__,
        "github_run_id": os.environ.get("GITHUB_RUN_ID"),
        "github_sha": os.environ.get("GITHUB_SHA"), "checks": {}, "passed": False,
        "blas_thread_limit": 1,
    }
    try:
        with tempfile.TemporaryDirectory(prefix="motifsurvey-check-") as temporary:
            copy = Path(temporary) / "source"
            shutil.copytree(root, copy, ignore=shutil.ignore_patterns(".git", "build", "dist", ".venv", "*.egg-info", "__pycache__", "data-cache", "check-results"))
            if args.cache:
                dataset = args.cache / "email-Eu-core.txt.gz"
                if dataset.exists():
                    shutil.copy2(dataset, copy / "validation/followup" / dataset.name)
            def execute(name, command):
                completed = subprocess.run([sys.executable, *command], cwd=copy, capture_output=True, text=True, encoding="utf-8", errors="replace", env=dict(os.environ, PYTHONUTF8="1"))
                log = completed.stdout + "\n" + completed.stderr
                (args.output / (name + ".log")).write_text(log, encoding="utf-8")
                report["checks"][name] = {"returncode": completed.returncode, "passed": completed.returncode == 0}
                if completed.returncode:
                    raise RuntimeError(f"{name} failed; see its log")
                return log
            units = execute("unit_tests", ["-m", "unittest", "discover", "-s", "tests", "-v"])
            report["checks"]["unit_tests"]["count"] = int(re.search(r"Ran (\d+) tests", units).group(1))
            cli = execute("cli", ["-m", "motifsurvey", "examples/tiny.json"])
            if json.loads(cli)["version"] != motifsurvey.__version__:
                raise RuntimeError("CLI version mismatch")
            execute("stress", ["validation/followup/stress.py"])
            stress = json.loads((copy / "validation/followup/stress_results.json").read_text())
            if stress["failures"] or any(c["result"] == "unexpected success" for c in stress["resource_checks"]):
                raise RuntimeError("Stress validation recorded a failure")
            report["checks"]["stress"].update(cases=len(stress["cases"]), expected_rejections=len(stress["resource_checks"]))
            execute("exact_and_replay", ["validation/followup/validate_v020.py"])
            validation = json.loads((copy / "validation/followup/v020_validation.json").read_text())
            if not validation["all_passed"]:
                raise RuntimeError("Exact/replay validation failed")
            report["checks"]["exact_and_replay"].update(exact_samples=validation["exact_samples"], public_replays=validation["public_replay_checks"], cached_native_comparisons=len(validation["native_fixture_checks"]))
            execute("dependencies", ["-m", "pip", "check"])
        report["passed"] = True
    except Exception as exc:
        report["error"] = str(exc)
        raise
    finally:
        (args.output / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
