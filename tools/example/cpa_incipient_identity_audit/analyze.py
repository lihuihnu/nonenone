"""Analyze diagnostic CSVs; completion is not a WLV-WL acceptance gate.

Usage: python analyze.py COARSE_OUTPUT HALF_STEP_OUTPUT ANALYSIS_OUTPUT
Requires numpy and pandas. No EOS parameters are read, fitted or changed.
"""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd


def analyze(coarse: Path, half: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    states = pd.read_csv(coarse / "incipient_identity.csv")
    a = pd.read_csv(coarse / "oil_tpd_hessian.csv")
    b = pd.read_csv(half / "oil_tpd_hessian.csv")
    keys = ["T_K", "offset_MPa"]
    if states.empty or states[keys].duplicated().any():
        raise ValueError("Empty or duplicated state keys")
    checks = []
    for key, group in a.groupby(keys):
        other = b[(b.T_K == key[0]) & (b.offset_MPa == key[1])]
        h = group.pivot(index="row", columns="column", values="H_central").to_numpy()
        h2 = other.pivot(index="row", columns="column", values="H_central").to_numpy()
        if h.shape != (4, 4) or h2.shape != (4, 4):
            raise ValueError(f"Incomplete Hessian at {key}")
        if not np.isfinite(h).all() or not np.isfinite(h2).all():
            raise ValueError(f"Nonfinite Hessian at {key}")
        checks.append({
            "T_K": key[0], "offset_MPa": key[1],
            "min_eigenvalue_h": float(np.linalg.eigvalsh((h + h.T) / 2)[0]),
            "min_eigenvalue_h2": float(np.linalg.eigvalsh((h2 + h2.T) / 2)[0]),
            "max_H_difference": float(np.max(np.abs(h - h2))),
            "max_H_asymmetry": float(np.max(np.abs(h2 - h2.T))),
        })
    checks = pd.DataFrame(checks)
    merged = states.merge(checks, on=keys, validate="one_to_one", how="left")
    if merged.min_eigenvalue_h.isna().any():
        raise ValueError("Some state rows lack Hessian diagnostics")
    merged.to_csv(out / "identity_analysis.csv", index=False)
    checks.to_csv(out / "step_refinement.csv", index=False)
    near = merged[np.isclose(merged.offset_MPa, 0.020)]
    near.to_csv(out / "summary_at_minus_0p02MPa.csv", index=False)
    summary = {
        "diagnostic_states": len(merged),
        "max_fd_eigenvalue_change": float(np.max(np.abs(
            checks.min_eigenvalue_h - checks.min_eigenvalue_h2))),
        "max_hessian_asymmetry": float(checks.max_H_asymmetry.max()),
        "max_candidate_fugacity_residual": float(states.candidate_fugacity_residual.max()),
        "identity_status": "UNRESOLVED_HIGH_T_NEAR_OIL_COMPOSITION_SPLIT",
        "full_figure7_acceptance": "BLOCKED",
        "note": "Eigenvalue sign diagnoses local composition curvature; neither a Gas slot nor a unique root proves vapor identity. No density-ratio threshold is used to force acceptance."
    }
    (out / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("usage: analyze.py COARSE_OUTPUT HALF_STEP_OUTPUT ANALYSIS_OUTPUT")
    analyze(*(Path(p) for p in sys.argv[1:]))
