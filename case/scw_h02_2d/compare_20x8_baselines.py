#!/usr/bin/env python3
"""Compare a new strict 20x8 H02 result against certified run49/run86 baselines."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path

SCOPE="CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION"
TARGETS=(1.0,2.0)
MODES=("A","B","C")
LIMITS={
    "RF_H_abs":0.001,          # 0.1 percentage point
    "cum_Heavy_rel":0.002,     # 0.2 %
    "deltaP_rel":0.01,         # 1 %
    "mu_rel":0.005,            # 0.5 %
    "mechanism_RF_abs":0.001,  # 0.1 percentage point
}

def read_rows(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"empty CSV: {path}")
    return rows

def load(path):
    out={}
    for r in read_rows(path):
        key=(r["mode"],float(r["target_pvi"]))
        out[key]=r
    missing=[(m,t) for m in MODES for t in TARGETS if (m,t) not in out]
    if missing:
        raise ValueError(f"missing rows in {path}: {missing}")
    return out

def val(r,k):
    x=float(r[k])
    if not math.isfinite(x):
        raise ValueError(f"nonfinite {k}")
    return x

def rel(a,b):
    return abs(a-b)/max(abs(a),abs(b),1e-30)

def write_csv(path,rows):
    if not rows:
        return
    with Path(path).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def compare_reference(name,current,reference):
    rows=[]; gates=[]
    for t in TARGETS:
        for m in MODES:
            c=current[(m,t)]; r=reference[(m,t)]
            metrics={
                "RF_H_abs":abs(val(c,"RF_H")-val(r,"RF_H")),
                "cum_Heavy_rel":rel(val(c,"cum_Heavy_kg"),val(r,"cum_Heavy_kg")),
                "deltaP_rel":rel(val(c,"deltaP_MPa"),val(r,"deltaP_MPa")),
                "mu_rel":rel(val(c,"mu_Heavy_carrier_Pa_s"),val(r,"mu_Heavy_carrier_Pa_s")),
            }
            rows.append({
                "scope":SCOPE,"reference":name,"target_pvi":t,"mode":m,
                "RF_H_reference":val(r,"RF_H"),"RF_H_current":val(c,"RF_H"),
                "RF_H_abs_difference":metrics["RF_H_abs"],
                "cum_Heavy_reference_kg":val(r,"cum_Heavy_kg"),
                "cum_Heavy_current_kg":val(c,"cum_Heavy_kg"),
                "cum_Heavy_relative_difference":metrics["cum_Heavy_rel"],
                "deltaP_reference_MPa":val(r,"deltaP_MPa"),
                "deltaP_current_MPa":val(c,"deltaP_MPa"),
                "deltaP_relative_difference":metrics["deltaP_rel"],
                "mu_reference_Pa_s":val(r,"mu_Heavy_carrier_Pa_s"),
                "mu_current_Pa_s":val(c,"mu_Heavy_carrier_Pa_s"),
                "mu_relative_difference":metrics["mu_rel"],
            })
            for metric,x in metrics.items():
                gates.append({
                    "scope":SCOPE,"reference":name,"target_pvi":t,"mode":m,
                    "metric":metric,"value":x,"limit":LIMITS[metric],
                    "status":"PASS" if x<=LIMITS[metric] else "FAIL",
                })

    mechanism=[]
    for t in TARGETS:
        for label,hi,lo in (("B_minus_A","B","A"),("C_minus_B","C","B")):
            cur=val(current[(hi,t)],"RF_H")-val(current[(lo,t)],"RF_H")
            ref=val(reference[(hi,t)],"RF_H")-val(reference[(lo,t)],"RF_H")
            diff=abs(cur-ref)
            mechanism.append({
                "scope":SCOPE,"reference":name,"target_pvi":t,
                "mechanism":label,"reference_delta_RF_H":ref,
                "current_delta_RF_H":cur,"absolute_difference":diff,
                "limit":LIMITS["mechanism_RF_abs"],
                "status":"PASS" if diff<=LIMITS["mechanism_RF_abs"] else "FAIL",
            })
    return rows,gates,mechanism

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--current",required=True,type=Path)
    p.add_argument("--run49",required=True,type=Path)
    p.add_argument("--run86",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)

    current=load(a.current)
    all_rows=[]; all_gates=[]; all_mech=[]
    for name,path in (("run49",a.run49),("run86",a.run86)):
        rows,gates,mech=compare_reference(name,current,load(path))
        all_rows+=rows; all_gates+=gates; all_mech+=mech

    write_csv(a.out/"matched_pvi_baseline_comparison.csv",all_rows)
    write_csv(a.out/"baseline_gates.csv",all_gates)
    write_csv(a.out/"mechanism_increment_comparison.csv",all_mech)

    ok=all(x["status"]=="PASS" for x in all_gates+all_mech)
    status={
        "scope":SCOPE,
        "baseline_consistency":"PASS" if ok else "FAIL",
        "references":["run49","run86"],
        "limits":LIMITS,
    }
    (a.out/"baseline_status.json").write_text(
        json.dumps(status,indent=2)+"\n",encoding="utf-8")

    print(json.dumps(status,indent=2))
    for name in ("run49","run86"):
        subset=[x for x in all_gates if x["reference"]==name]
        worst=max(subset,key=lambda x:x["value"]/x["limit"])
        print("worst",name,worst)
    for x in all_mech:
        print("mechanism",x)

    raise SystemExit(0 if ok else 2)

if __name__=="__main__":
    main()
