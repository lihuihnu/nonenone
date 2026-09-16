# SCW–light/middle/heavy kerogen-product displacement

This is Stage 2 of the SCW/kerogen-product flow experiment.  It keeps the
Stage-1 grid, rock, temperature, pressure, wells, rate and PVI fixed, and
replaces the single heavy proxy with three product lumps:

| Lump | Numerical representative | Role |
|---|---|---|
| light | nC4 | volatile cracking product |
| middle | nC10 | mobile liquid/intermediate product |
| heavy | squalane (branched C30H62) | heavy residual/product proxy |

The current overall mole fractions are H2O/light/middle/heavy =
0.25/0.078125/0.1953125/0.4765625.  PR, SW and CPA all flash this state as one
hydrocarbon-rich phase with `Sw=0`; injection creates the water-rich phase.
At 653.2 K the H2O–nC4/nC10/squalane non-aqueous BIPs are respectively
0.5091/0.2618373654/0.0532336595.  Their evidence levels are published prior,
temperature extrapolation, and target-window composition fit.  Hydrocarbon
cross-BIPs remain zero, so this is still a mechanism-screening fluid rather
than a characterized kerogen assay.

```bash
make case CASE=scw_kerogen_lmh_1d -j
make run CASE=scw_kerogen_lmh_1d NP=1 \
  RESULT_DIR=results/base
python case/scw_kerogen_common/analyze_results.py \
  case/scw_kerogen_lmh_1d/results/base
```

Use `producer_composition.csv` to compare water breakthrough and
light/middle/heavy mass fractions.  Use `viscosity_history.csv` together with
phase saturations; a phase-role transition can move the hydrocarbon-enriched
fluid between the oil and gas slots at strict-SCW conditions.

The previously checked 1-PVI result below belongs to the superseded initial
composition and BIPs.  It reached 4 days with maximum whole-run
component-balance relative error `3.75e-8`.  The producer H2O mass fraction reaches about 0.860;
within the produced hydrocarbon fraction, light/middle/heavy change from about
0.0194/0.1189/0.8618 at the first solved output to
0.0409/0.1198/0.8393 at 1 PVI.  These numbers are historical and must not be
reported as results of the current baseline.

The full staged design and acceptance criteria are in
[`../../docs/SCW_KEROGEN_FLOW_EXPERIMENT.md`](../../docs/SCW_KEROGEN_FLOW_EXPERIMENT.md).
