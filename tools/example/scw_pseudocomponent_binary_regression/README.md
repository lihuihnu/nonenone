# Generic H2O–pseudo-component zero-flow PR regression harness

This executable is the common zero-flow regression/validation harness for the four experimental boiling-range pseudo-components in `case/scw_kerogen_lumped_flow_study`.

It reuses production `CubicEquationOfState` and `CubicThreePhaseFlash`; it does **not** implement a separate flash model.

## Inputs

```text
case/scw_kerogen_lumped_flow_study/binary_pr_calibration/binary_systems.csv
case/scw_kerogen_lumped_flow_study/binary_pr_calibration/binary_observations.csv
```

`binary_systems.csv` supplies the current pseudo-component `Tc/Pc/Vc/omega/MW`, the `kij(T)` reference temperature and parameter bounds, plus a per-system fit-enable policy.

`binary_observations.csv` uses one common schema for coexistence, boundary and density/volume validation rows. Missing experimental quantities are explicit `nan`; they are never silently replaced by model values.

## BIP regression

The only parameters fitted by this harness are

```text
kij(T) = k_ref + b * (1/T - 1/T_ref)
```

and `b` is fitted only when calibration rows span at least two temperatures. With only one calibration temperature the harness performs a constant-`kij` diagnostic.

The BIP objective uses **only both-branch coexistence compositions** from rows labelled `calibration`. Density or molar-volume targets are never put into the BIP objective; this prevents `kij` from compensating for a density-model error.

## Independent gates

For every system the output `gate_summary.csv` separately evaluates:

- phase count;
- both coexistence compositions;
- phase density/volume targets;
- boundary rows;
- held-out validation rows.

A system is `PASS` only if every gate has at least one experimental row and all rows pass. Missing data therefore produces `NO_DATA` and keeps the reservoir gate blocked.

`OIL_HEAVY` has `fit_enabled=0` until an admissible heavy-oil proxy envelope is available near the target pressure. This is intentional: Heavy may not be tuned merely to make flash converge.

## Current data state

The exact Stevenson water–squalane tie-lines already in the repository seed the `OIL_MIDDLE` proxy diagnostic. n-hexane and dodecane are still data-blocked because publicly recovered sources currently confirm dataset membership/counts/ranges but not the original pointwise tables. Recovery status is documented in:

```text
case/scw_kerogen_lumped_flow_study/binary_pr_calibration/raw_recovery/
```

## Run

From the repository root:

```bash
make -C tools run-scw-pseudocomponent-binary-regression CXX=clang++
```

The target refuses to overwrite its existing result directory. A nonzero exit code `3` means the executable ran correctly but the reservoir-entry gate remains blocked; it is not a numerical-crash code.
