#!/usr/bin/env python3
"""Compare strict H02 20x8 dt and 60x20 grid convergence at matched actual PVI.

All inputs must already pass analyze_h02.py independently. This script adds
pre-registered numerical screening gates; it does not change physical models or
claim experimental validation.
"""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path

SCOPE="CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION"
TARGETS=(1.0,2.0)
MODES=("A","B","C")
TIME_LIMITS=dict(rf_abs=0.001, dp_rel=0.01, heavy_rel=0.002, mu_rel=0.005)
GRID_LIMITS=dict(rf_abs=0.005, dp_rel=0.05, heavy_rel=0.01, mu_rel=0.01)

def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows: raise ValueError(f"empty table: {path}")
    return rows

def num(row,key):
    value=float(row[key])
    if not math.isfinite(value): raise ValueError(f"nonfinite {key}")
    return value

def load_summary(root):
    rows=read_csv(Path(root)/"analysis"/"h02_mechanism_at_1_2_pvi.csv")
    out={}
    for r in rows:
        key=(r["mode"],float(r["target_pvi"]))
        out[key]=r
    missing=[(m,t) for t in TARGETS for m in MODES if (m,t) not in out]
    if missing: raise ValueError(f"missing matched-PVI rows: {missing}")
    return out

def rel_diff(a,b):
    scale=max(abs(a),abs(b),1e-30)
    return abs(b-a)/scale

def compare(label,base,test,limits):
    rows=[]; gates=[]
    for target in TARGETS:
        for mode in MODES:
            a,b=base[(mode,target)],test[(mode,target)]
            rf0,rf1=num(a,"RF_H"),num(b,"RF_H")
            dp0,dp1=num(a,"deltaP_MPa"),num(b,"deltaP_MPa")
            h0,h1=num(a,"cum_Heavy_kg"),num(b,"cum_Heavy_kg")
            mu0,mu1=num(a,"mu_Heavy_carrier_Pa_s"),num(b,"mu_Heavy_carrier_Pa_s")
            w0,w1=num(a,"cum_Heavy_Water_kg"),num(b,"cum_Heavy_Water_kg")
            wf0,wf1=num(a,"water_carried_fraction"),num(b,"water_carried_fraction")
            rec=dict(scope=SCOPE,comparison=label,mode=mode,target_pvi=target,
                     RF_H_base=rf0,RF_H_test=rf1,RF_H_abs_difference=abs(rf1-rf0),
                     deltaP_base_MPa=dp0,deltaP_test_MPa=dp1,deltaP_relative_difference=rel_diff(dp0,dp1),
                     cum_Heavy_base_kg=h0,cum_Heavy_test_kg=h1,cum_Heavy_relative_difference=rel_diff(h0,h1),
                     mu_base_Pa_s=mu0,mu_test_Pa_s=mu1,mu_relative_difference=rel_diff(mu0,mu1),
                     cum_Heavy_Water_base_kg=w0,cum_Heavy_Water_test_kg=w1,
                     water_carried_fraction_base=wf0,water_carried_fraction_test=wf1)
            rows.append(rec)
            checks={
              "RF_H_abs":(rec["RF_H_abs_difference"],limits["rf_abs"]),
              "deltaP_rel":(rec["deltaP_relative_difference"],limits["dp_rel"]),
              "cum_Heavy_rel":(rec["cum_Heavy_relative_difference"],limits["heavy_rel"]),
              "mu_rel":(rec["mu_relative_difference"],limits["mu_rel"]),
            }
            for metric,(value,limit) in checks.items():
                gates.append(dict(scope=SCOPE,comparison=label,mode=mode,target_pvi=target,
                                  metric=metric,value=value,limit=limit,
                                  status="PASS" if value<=limit else "FAIL"))
    # Mechanism increments are reported separately; no extra threshold is invented.
    mech=[]
    for target in TARGETS:
        for name,x,y in (("B_minus_A","B","A"),("C_minus_B","C","B")):
            for tag,data in (("base",base),("test",test)):
                rx=num(data[(x,target)],"RF_H")-num(data[(y,target)],"RF_H")
                dp=num(data[(x,target)],"deltaP_MPa")-num(data[(y,target)],"deltaP_MPa")
                mu=num(data[(x,target)],"mu_Heavy_carrier_Pa_s")-num(data[(y,target)],"mu_Heavy_carrier_Pa_s")
                mech.append(dict(scope=SCOPE,comparison=label,target_pvi=target,
                                 mechanism=name,dataset=tag,delta_RF_H=rx,
                                 delta_deltaP_MPa=dp,delta_mu_Pa_s=mu))
    return rows,gates,mech

def write_csv(path,rows):
    if not rows:return
    with Path(path).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--dt2",required=True,type=Path,help="20x8, dt_max=2s ensemble root")
    p.add_argument("--dt1",required=True,type=Path,help="20x8, dt_max=1s ensemble root")
    p.add_argument("--grid60",required=True,type=Path,help="60x20, dt_max=2s ensemble root")
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    dt2=load_summary(a.dt2);dt1=load_summary(a.dt1);fine=load_summary(a.grid60)
    tr,tg,tm=compare("20x8_dt2_vs_dt1",dt2,dt1,TIME_LIMITS)
    gr,gg,gm=compare("20x8_vs_60x20_dt2",dt2,fine,GRID_LIMITS)
    gates=tg+gg
    write_csv(a.out/"time_step_convergence.csv",tr)
    write_csv(a.out/"grid_convergence.csv",gr)
    write_csv(a.out/"convergence_gates.csv",gates)
    write_csv(a.out/"mechanism_increment_stability.csv",tm+gm)
    status={
      "scope":SCOPE,
      "time_step_gate":"PASS" if all(x["status"]=="PASS" for x in tg) else "FAIL",
      "grid_gate":"PASS" if all(x["status"]=="PASS" for x in gg) else "FAIL",
      "time_limits":TIME_LIMITS,
      "grid_limits":GRID_LIMITS,
      "water_carriage_gate":"REPORT_ONLY_TRACE_SCALE",
      "physical_validation":"NOT_VALIDATED",
    }
    (a.out/"convergence_status.json").write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(status,indent=2))
    print("\nWorst failed/passing margins:")
    for label,subset in (("TIME",tg),("GRID",gg)):
        worst=max(subset,key=lambda x:x["value"]/x["limit"] if x["limit"] else math.inf)
        print(label,worst)
    raise SystemExit(0 if status["time_step_gate"]=="PASS" and status["grid_gate"]=="PASS" else 2)

if __name__=="__main__":
    main()
