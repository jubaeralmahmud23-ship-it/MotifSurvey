"""Verify the actual native R research-workflow runs; MIT License."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def check(folder):
    native = json.loads((folder / "native_results.json").read_text())
    expected = {x["name"]: x for x in json.loads((folder / "native_expected.json").read_text())}
    rows = []
    assert len(native["results"]) == 3
    for result in native["results"]:
        target = expected[result["name"]]
        totals = np.asarray(result["totals"])
        covariance = np.asarray(result["covariance"])
        np.testing.assert_allclose(totals, target["result"]["estimated_totals"], rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(covariance, target["result"]["estimated_covariance"], rtol=1e-8, atol=1e-8)
        errors = []
        for profile in target["profile"]["profile"]:
            c = np.array([0, -profile["threshold"], 3])
            actual = np.array([c @ totals, c @ covariance @ c])
            wanted = np.array([profile["estimate"], profile["estimated_variance"]])
            np.testing.assert_allclose(actual, wanted, rtol=1e-8, atol=1e-8)
            errors.extend((abs(actual-wanted)/np.maximum(1, abs(wanted))).tolist())
        rows.append({"allocation": result["name"], "motif_units": result["units"],
                     "subset_keys": target["result"]["diagnostics"]["subset_keys"],
                     "maximum_scaled_profile_error": max(errors), "passed": True})
    report = {"scope": "Three actual R survey runs on the first evaluation sample under each allocation; nine contrast/variance pairs. Not independent datasets or a new timing benchmark.",
              "R": native["R"], "versions": native["versions"], "checks": rows,
              "native_adapter_sha256": hashlib.sha256((Path(__file__).resolve().parents[1] / "followup/native_compare.R").read_bytes()).hexdigest(),
              "maximum_scaled_profile_error": max(x["maximum_scaled_profile_error"] for x in rows), "passed": True}
    (folder / "native_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    check(parser.parse_args().folder)
