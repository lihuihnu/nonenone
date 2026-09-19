#!/usr/bin/env python3
import csv, math, sys
from pathlib import Path

def rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))

def f(x):
    try: return float(x)
    except: return float("nan")

root=Path(sys.argv[1])
out=Path(sys.argv[2])
out.parent.mkdir(parents=True, exist_ok=True)
summary=[]
for case_dir in sorted(p for p in root.iterdir() if p.is_dir()):
    mb=rows(case_dir/"component_mass_balance.csv")
    nl=rows(case_dir/"nonlinear_solve_history.csv")
    ts=rows(case_dir/"time_step_history.csv")
    step1=[r for r in mb if int(r["step"])==1]
    errs={r["component"]: f(r["balance_error_kg"]) for r in step1}
    inv={r["component"]: f(r["initial_inventory_kg"]) for r in step1}
    inj={r["component"]: f(r["cumulative_injected_kg"]) for r in step1}
    scaled={}
    for comp,e in errs.items():
        denom=inv.get(comp,0.0)+inj.get(comp,0.0)
        scaled[comp]=abs(e)/denom if denom>1e-30 else float("nan")
    maxres=max(f(r["final_residual"]) for r in nl)
    worst=max(nl,key=lambda r:f(r["final_residual"]))
    rejected=sum(1 for r in ts if r["status"]=="REJECT")
    summary.append({
      "case":case_dir.name,
      "heavy_error_kg":errs.get("OIL_HEAVY",float("nan")),
      "heavy_relative_error":scaled.get("OIL_HEAVY",float("nan")),
      "water_error_kg":errs.get("H2O",float("nan")),
      "water_relative_error":scaled.get("H2O",float("nan")),
      "max_final_residual":maxres,
      "worst_solve_id":worst["solve_id"],
      "worst_reason":worst["reason"],
      "accepted_steps":sum(1 for r in ts if r["status"]=="ACCEPT"),
      "rejected_steps":rejected,
    })
with out.open("w",newline="") as fobj:
    w=csv.DictWriter(fobj,fieldnames=list(summary[0].keys()))
    w.writeheader(); w.writerows(summary)
print(out.read_text())
