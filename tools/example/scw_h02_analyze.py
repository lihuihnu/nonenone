#!/usr/bin/env python3
"""Summarize SCW-H02 A/B/C at matched actual PVI."""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import pandas as pd
import numpy as np

def rf_column(df):
    for name in ("RF_OIL_HEAVY", "RF_total_hydrocarbon"):
        if name in df.columns:
            return name
    candidates=[c for c in df.columns if c.startswith("RF_") and "H2O" not in c]
    if not candidates:
        raise RuntimeError("No Heavy recovery column found")
    return candidates[0]

def interp(df, col, pvi):
    x=df["pvi"].to_numpy(float)
    y=df[col].to_numpy(float)
    mask=np.isfinite(x)&np.isfinite(y)
    x=x[mask]; y=y[mask]
    if len(x)<2 or pvi<x.min()-1e-12 or pvi>x.max()+1e-12:
        return math.nan
    return float(np.interp(pvi,x,y))

def read_mode(root, mode, target):
    d=root/mode
    prod=pd.read_csv(d/"producer_composition.csv")
    mb=pd.read_csv(d/"component_mass_balance.csv")
    wh=pd.read_csv(d/"well_history.csv")
    rd=pd.read_csv(d/"reservoir_diagnostics.csv")
    ts=pd.read_csv(d/"time_step_history.csv")
    rf=rf_column(prod)
    last=prod.iloc[-1]
    inj=wh[wh["type"]=="INJECTOR"]
    result={
      "mode":mode,
      "final_time_s":float(last["time_s"]),
      "final_pvi":float(last["pvi"]),
      "target_reached":bool(float(last["pvi"])>=target-1e-8),
      "rf_heavy_final":float(last[rf]),
      "rf_heavy_pvi_1":interp(prod,rf,1.0),
      "rf_heavy_pvi_2":interp(prod,rf,2.0),
      "max_abs_component_relative_error":float(mb["relative_error"].abs().max()),
      "max_injector_bhp_bar":float(inj["bhp"].max()),
      "final_injector_control":str(inj.iloc[-1]["active_control"]),
      "accepted_internal_steps":int((ts["status"]=="ACCEPT").sum()),
      "rejected_internal_steps":int((ts["status"]=="REJECT").sum()),
      "pressure_min_bar":float(rd.iloc[-1]["pressure_min"]),
      "pressure_avg_bar":float(rd.iloc[-1]["pressure_avg"]),
      "pressure_max_bar":float(rd.iloc[-1]["pressure_max"]),
      "oil_viscosity_min_Pa_s":float(rd.iloc[-1]["mu_o_min"]),
      "oil_viscosity_avg_Pa_s":float(rd.iloc[-1]["mu_o_avg"]),
      "water_saturation_max":float(rd.iloc[-1]["s_w_max"]),
    }
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("root",type=Path)
    ap.add_argument("--target-pvi",type=float,default=2.0)
    ap.add_argument("--require-target",action="store_true")
    ap.add_argument("--max-balance-error",type=float,default=1e-5)
    args=ap.parse_args()
    rows=[read_mode(args.root,m,args.target_pvi) for m in "ABC"]
    by={r["mode"]:r for r in rows}
    comparisons={}
    for p in (1,2):
        a=by["A"][f"rf_heavy_pvi_{p}"]
        b=by["B"][f"rf_heavy_pvi_{p}"]
        c=by["C"][f"rf_heavy_pvi_{p}"]
        comparisons[f"pvi_{p}"]={
          "RF_A":a,"RF_B":b,"RF_C":c,
          "delta_exchange_B_minus_A":b-a if math.isfinite(a) and math.isfinite(b) else math.nan,
          "delta_viscosity_C_minus_B":c-b if math.isfinite(b) and math.isfinite(c) else math.nan,
        }
    out={"modes":rows,"comparisons":comparisons}
    args.root.mkdir(parents=True,exist_ok=True)
    (args.root/"h02_summary.json").write_text(json.dumps(out,indent=2,allow_nan=True)+"\n")
    pd.DataFrame(rows).to_csv(args.root/"h02_summary.csv",index=False)
    print(json.dumps(out,indent=2,allow_nan=True))
    bad=[r for r in rows if r["max_abs_component_relative_error"]>args.max_balance_error]
    if bad:
        raise SystemExit("mass-balance certificate failed: "+",".join(r["mode"] for r in bad))
    if any(r["max_injector_bhp_bar"]>300.0001 for r in rows):
        raise SystemExit("injector BHP exceeded 30 MPa design cap")
    if args.require_target and any(not r["target_reached"] for r in rows):
        raise SystemExit("one or more modes did not reach target actual PVI")
if __name__=="__main__":
    main()
