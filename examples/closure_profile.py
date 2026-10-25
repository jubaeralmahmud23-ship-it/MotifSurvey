"""Reuse one MotifSurvey result for prespecified closure thresholds.

MIT License. This is linear contrast postprocessing, not a test or interval.
Thresholds selected after inspecting results must be described as exploratory.
"""
import argparse
import json
from numbers import Real
from decimal import Decimal
from pathlib import Path

import numpy as np


def closure_profile(result, thresholds):
    """Return estimates and raw variances without re-enumerating motifs."""
    if result.get("motif_order") != ["edges", "wedges", "triangles"]:
        raise ValueError("Expected MotifSurvey edge, wedge, triangle ordering")
    totals = np.asarray(result["estimated_totals"], dtype=float)
    covariance = np.asarray(result["estimated_covariance"], dtype=float)
    if totals.shape != (3,) or covariance.shape != (3, 3):
        raise ValueError("Expected three totals and a 3 by 3 covariance matrix")
    if not np.isfinite(totals).all() or not np.isfinite(covariance).all():
        raise ValueError("Totals and covariance must be finite")
    if not np.allclose(covariance, covariance.T, rtol=1e-12, atol=1e-12):
        raise ValueError("Covariance must be symmetric")
    values = list(thresholds)
    if not values or any(isinstance(t, (bool, np.bool_)) or not isinstance(t, (Real, Decimal)) for t in values):
        raise ValueError("Supply numeric thresholds between zero and one")
    try:
        values = np.asarray(values, dtype=float)
    except (ValueError, TypeError, OverflowError):
        raise ValueError("Supply numeric thresholds between zero and one") from None
    if values.ndim != 1 or not np.isfinite(values).all() or np.any((values < 0) | (values > 1)):
        raise ValueError("Supply numeric thresholds between zero and one")
    rows = []
    for threshold in values:
        c = np.array([0.0, -threshold, 3.0])
        variance = float(c @ covariance @ c)
        rows.append({
            "threshold": float(threshold), "estimate": float(c @ totals),
            "estimated_variance": variance,
            "diagonal_only_variance": float((c * c) @ np.diag(covariance)),
            "covariance_contribution": float(-6 * threshold * covariance[1, 2]),
            "negative_estimated_variance": variance < 0,
        })
    return {
        "input_sha256": result.get("input_sha256"),
        "software_version": result.get("version"),
        "source_warnings": list(result.get("warnings", [])),
        "interpretation": "D(t) = 3*T_triangle - t*T_wedge. Raw estimated variances are retained. No significance, confidence interval, simultaneous coverage or threshold-selection guarantee.",
        "profile": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--thresholds", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = closure_profile(json.loads(args.input.read_text(encoding="utf-8")), args.thresholds)
    content = json.dumps(output, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(content)
    else:
        print(content, end="")
