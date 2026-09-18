#!/usr/bin/env python3
"""Audit H02 native output. Conditional study, not real-Heavy validation.

Uses accepted-step cumulative ledgers, not reintegration of sparse output rates.
No extrapolation, no missing-column zero fill, no 1-kg relative-error floor.
Viscosity is reconstructed from the existing H02 assumed law at one producer
perforation, not an experimentally validated mixture viscosity measurement.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

PHASES = ("Oil", "Gas", "Water")
SCOPE = "CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION"
MASS_LIMIT = 1e-6

def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"Empty required table: {path}")
    return rows

def num(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite required value: {key}={value}")
    return value

def write_csv(path, rows):
    if rows:
        with Path(path).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

def bracket(rows, target):
    if not rows or not math.isfinite(target):
        raise ValueError("Invalid target or empty series")
    x = [num(row, "pvi") for row in rows]
    if any(b <= a for a, b in zip(x, x[1:])):
        raise ValueError("PVI must be strictly increasing")
    if target < x[0] or target > x[-1]:
        raise ValueError(f"Target {target} outside actual coverage [{x[0]}, {x[-1]}]")
    for i, value in enumerate(x):
        if target == value:
            return rows[i], rows[i], 0.0
        if value > target:
            return rows[i-1], rows[i], (target-x[i-1])/(value-x[i-1])
    raise AssertionError("Unreachable bracket")

def interpolate(rows, target, key):
    a, b, fraction = bracket(rows, target)
    va, vb = float(a[key]), float(b[key])
    if fraction == 0.0:
        return va
    # Do not bridge undefined/no-flow intervals or fill outside coverage.
    return va + fraction*(vb-va) if math.isfinite(va) and math.isfinite(vb) else math.nan

def h02_mu(w, mode):
    if not math.isfinite(w) or not 0.0 <= w <= 1.0:
        raise ValueError("Invalid water mass fraction")
    if mode == "B":
        s = min(1.0, max(0.0, (w-0.25)/0.5))
        v = s*s*(3.0-2.0*s)
    elif mode == "C":
        v = w
    else:
        raise ValueError("Composition law is only defined for B/C")
    return math.exp(math.log(0.002)+v*math.log(5e-5/0.002))

def load_mode(root, mode):
    raw = read_csv(root/mode/"producer_composition.csv")
    wells = read_csv(root/mode/"well_history.csv")
    lookup = {(int(row["step"]), row["name"]): row for row in wells}
    rows, rates_error, cumulative_error = [], 0.0, 0.0
    for row in raw:
        step = int(row["step"])
        inj, prod = lookup[(step,"SCW_INJ")], lookup[(step,"PROD")]
        if int(prod["perforations"]) != 1:
            raise ValueError("Viscosity reconstruction requires one producer perforation")
        if num(prod,"m_component_OIL_HEAVY") > 1e-14:
            raise ValueError("Producer reverse flow requires gross phase ledgers")
        rec = dict(scope=SCOPE, mode=mode, step=step,
            time_s=num(row,"time_s"), pvi=num(row,"pvi"),
            RF_H=num(row,"RF_OIL_HEAVY"),
            cum_Heavy_kg=num(row,"cum_OIL_HEAVY_produced_kg"),
            deltaP_MPa=(num(inj,"bhp")-num(prod,"bhp"))*0.1,
            injector_BHP_MPa=num(inj,"bhp")*0.1,
            producer_BHP_MPa=num(prod,"bhp")*0.1,
            injector_control=inj["active_control"])
        weighted, total_rate, phase_cum = 0.0, 0.0, 0.0
        for phase in PHASES:
            rate = num(row,f"m_dot_{phase}_OIL_HEAVY_produced_kg_s")
            cum = num(row,f"cum_{phase}_OIL_HEAVY_produced_kg")
            if rate < 0 or cum < 0:
                raise ValueError("Negative production magnitude")
            if mode == "A":
                mu = 5e-5 if phase == "Water" else 0.002
            else:
                water = num(row,f"m_dot_{phase}_H2O_produced_kg_s")
                mu = h02_mu(water/(water+rate), mode) if water+rate > 0 else math.nan
            rec[f"m_dot_Heavy_{phase}_kg_s"] = rate
            rec[f"cum_Heavy_{phase}_kg"] = cum
            rec[f"mu_{phase}_Pa_s"] = mu
            if rate > 0:
                weighted += rate*mu
            total_rate += rate
            phase_cum += cum
        rec["mu_Heavy_carrier_Pa_s"] = weighted/total_rate if total_rate > 0 else math.nan
        rates_error = max(rates_error,abs(total_rate-num(row,"m_dot_OIL_HEAVY_produced_kg_s")))
        cumulative_error = max(cumulative_error,abs(phase_cum-rec["cum_Heavy_kg"]))
        rows.append(rec)
    bracket(rows,1.0)
    bracket(rows,2.0)
    return rows, rates_error, cumulative_error

def audit(root, out, mode_names="ABC"):
    out.mkdir(parents=True,exist_ok=True)
    requested=tuple(mode_names)
    if not requested or any(m not in "ABC" for m in requested):
        raise ValueError("modes must be a non-empty subset of A/B/C")
    if len(set(requested)) != len(requested):
        raise ValueError("modes must not contain duplicates")
    modes, gates, mass_summary = {}, [], []
    def gate(mode,name,value,limit):
        passed = math.isfinite(value) and value <= limit
        gates.append(dict(mode=mode,gate=name,value=value,limit=limit,status="PASS" if passed else "FAIL"))
    for mode in requested:
        rows, rate_err, cumulative_err = load_mode(root,mode)
        modes[mode] = rows
        gate(mode,"phase_component_rate_closure_kg_s",rate_err,1e-14)
        gate(mode,"phase_component_cumulative_closure_kg",cumulative_err,1e-12)
        if mode == "A":
            gate(mode,"no_transfer_water_carriage_kg",max(r["cum_Heavy_Water_kg"] for r in rows),1e-12)
        mb = read_csv(root/mode/"component_mass_balance.csv")
        components = {row["component"] for row in mb}
        if components != {"H2O","OIL_HEAVY"}:
            raise ValueError("Both conserved components must be independently audited")
        for comp in sorted(components):
            selected = [r for r in mb if r["component"] == comp]
            values = []
            for row in selected:
                initial = num(row,"initial_inventory_kg")
                injected = num(row,"cumulative_injected_kg")
                produced = num(row,"cumulative_produced_kg")
                current = num(row,"current_inventory_kg")
                error = current-initial-injected+produced
                budget = initial+injected
                # A zero-inventory component is scaled by actual injected mass;
                # a negligible empty state is handled explicitly, never by 1 kg.
                rel = abs(error)/budget if budget > 1e-24 else (0.0 if abs(error)<1e-24 else math.inf)
                values.append((rel,abs(error),int(row["step"])))
            worst = max(values,key=lambda v:v[0])
            gate(mode,f"{comp}_mass_relative",worst[0],MASS_LIMIT)
            mass_summary.append(dict(mode=mode,component=comp,
                max_relative_error=worst[0],worst_relative_step=worst[2],
                max_absolute_error_kg=max(v[1] for v in values),
                normalization="abs(M-M0-Min+Mout)/(M0+Min)"))
    summary = []
    for target in (1.0,2.0):
        for mode in requested:
            rows = modes[mode]
            a,b,f = bracket(rows,target)
            rec = dict(scope=SCOPE,target_pvi=target,mode=mode,
                       lower_actual_pvi=a["pvi"],upper_actual_pvi=b["pvi"],fraction=f,
                       method="BRACKETED_LINEAR_NO_EXTRAPOLATION")
            for key in ("time_s","RF_H","cum_Heavy_kg","cum_Heavy_Oil_kg","cum_Heavy_Gas_kg",
                        "cum_Heavy_Water_kg","deltaP_MPa","mu_Heavy_carrier_Pa_s"):
                rec[key] = interpolate(rows,target,key)
            rec["water_carried_fraction"] = rec["cum_Heavy_Water_kg"]/rec["cum_Heavy_kg"]
            summary.append(rec)
    bc = []
    if "B" in modes and "C" in modes:
        for i in range(101):
            pvi = i/50
            b = interpolate(modes["B"],pvi,"mu_Heavy_carrier_Pa_s")
            c = interpolate(modes["C"],pvi,"mu_Heavy_carrier_Pa_s")
            bc.append(dict(scope=SCOPE,pvi=pvi,mu_B_Heavy_carrier_Pa_s=b,
                           mu_C_Heavy_carrier_Pa_s=c,delta_mu_C_minus_B_Pa_s=c-b,
                           relative_delta_C_minus_B=(c-b)/b if b>0 else math.nan))
    write_csv(out/"h02_mechanism_timeseries.csv",[r for m in requested for r in modes[m]])
    write_csv(out/"h02_mechanism_at_1_2_pvi.csv",summary)
    write_csv(out/"h02_bc_viscosity_delta.csv",bc)
    write_csv(out/"h02_mass_audit.csv",mass_summary)
    write_csv(out/"h02_numerical_gates.csv",gates)
    status = dict(scope=SCOPE,pvi_coverage="PASS",mass_limit=MASS_LIMIT,
                  numerical_gate="PASS" if all(g["status"]=="PASS" for g in gates) else "FAIL",
                  physical_validation="NOT_VALIDATED",grid_time_convergence="NOT_ESTABLISHED")
    (out/"audit_status.json").write_text(json.dumps(status,indent=2)+"\n")
    return status

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    parser.add_argument("--modes",default="ABC")
    args = parser.parse_args()
    try:
        status = audit(args.root,args.out,args.modes)
    except (OSError,KeyError,ValueError) as exc:
        print(f"H02_DATA_INTEGRITY=FAIL: {exc}")
        raise SystemExit(2)
    print(json.dumps(status,indent=2))
    raise SystemExit(0 if status["numerical_gate"]=="PASS" else 2)

if __name__ == "__main__":
    main()
