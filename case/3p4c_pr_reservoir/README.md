# 3p4c_pr_reservoir

This is the first reservoir case that exercises the project's
`FullyCompositionalThreePhase` path end-to-end.

## Fluid and initialization

Components are ordered as:

1. H2O
2. CO2
3. CH4
4. nC16 (heavy hydrocarbon representative)

The user-supplied reservoir initial state is only:

- pressure: 200 bar
- temperature: 350 K
- overall mole fraction: `z = [0.65, 0.30, 0.02, 0.03]`

No phase saturation and no phase composition is entered in the case file.
`NaturalPetscRuntime::initializeUniformFromPTZ()` performs the PR flash and
writes the resulting O/G/W state into the primary solution and phase-state
vectors.

The standalone case regression (`three_phase_reservoir_case_test`) gives the
following expected initial flash (small differences at the last digits are
acceptable):

| Quantity | Oil-rich | Gas | Water-rich |
|---|---:|---:|---:|
| phase mole fraction beta | 0.16509154 | 0.18705867 | 0.64784979 |
| pore-volume saturation S | 0.37787066 | 0.31742857 | 0.30470078 |

Phase mole fractions `[H2O, CO2, CH4, nC16]` are approximately:

- oil-rich: `[0.00572568, 0.78312272, 0.05370433, 0.15744726]`
- gas: `[0.00969491, 0.90936764, 0.05951750, 0.02141996]`
- water-rich: `[0.99906064, 0.00093843, 9.324e-7, ~0]`

This demonstrates mutual partitioning rather than an immiscible-water/Henry
shortcut: H2O is finite in the oil/gas phases and CO2 is finite in the
water-rich phase.  The complete expected state (including Z and molar density)
is stored in `reference_initial_flash.csv` for quick comparison with cluster
`solution_step_0.csv` / `phase_state_step_0.csv` output.

## Reservoir problem

- structured grid: 30 x 15 x 1
- physical size: 1200 x 600 x 12 m
- porosity: 0.20
- permeability: 100/100/10 mD
- injector: left-center, pure CO2, total-rate control, 100 m3/day
- producer: right-center, BHP = 195 bar
- initial pressure: 200 bar
- default test horizon: 20 x 0.25 day = 5 days

This is intentionally a mild smoke problem.  Its first purpose is to test the
coupled reservoir path (component conservation + three Darcy phase fluxes +
three-phase fugacity equations + wells + phase switching) before increasing
rates, time horizon or geological complexity.

## Local thermodynamic/physics preflight

This check does not require PETSc:

```bash
cd ../../test
make bin/three_phase_reservoir_case_test CXX=g++
./bin/three_phase_reservoir_case_test
```

It verifies the exact case configuration through P-T-z flash, phase identity,
material/fugacity closure, EOS density/viscosity/mobility, component
accumulation, a pressure-driven face flux, and injector/producer source signs.
It also contains the regression for the water-rich trace-nC16 boundary: the
`N-1` primary-composition reconstruction must not create a negative dependent
mole fraction or a spurious `1e29+` oil-water fugacity residual.

## Build and run on the PETSc cluster

From `case/`:

```bash
make 3p4c_pr_reservoir -j
make reset-run CASE=3p4c_pr_reservoir
cd 3p4c_pr_reservoir
./run.sh
```

For a very short first integration test, override the case defaults without
editing C++:

```bash
./bin/3p4c_pr_reservoir -numSteps 2 -dt 0.05 -out_step 1
```

All standard result files are CSV.

## Scope of the fluid data

The case uses real component identities and a PR interaction pattern adapted
from published CO2/water/hydrocarbon three-phase benchmark work.  It is a
thermodynamic/numerical benchmark and is not claimed to be a field-fluid PVT
fit.  Before field prediction, Tc/Pc/omega/pseudo-component properties and BIPs
should be regressed to the target fluid's PVT/phase-equilibrium data.

## Per-component global conservation audit

From v22 this case enables `component_mass_balance.csv`. The ledger uses every
accepted internal backward-Euler step, not just `OUT_STEP` states. For each of
H2O/CO2/CH4/nC16 it reports initial/current reservoir inventory, cumulative
injection, cumulative production, expected inventory, absolute error and
relative error. Rejected adaptive attempts are excluded. The ledger definition is
documented in `../../output/include/output/metrics/MODULE.md`.
