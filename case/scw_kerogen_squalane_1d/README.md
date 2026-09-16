# SCW–single-heavy kerogen-product displacement

This is Stage 1 of the SCW/kerogen-product flow experiment.  It converts the
existing H2O–squalane strict-SCW flash fluid into a 60-cell isothermal
displacement case.  The reservoir is initialized from P-T-z; no phase
saturation or phase composition is prescribed.  Pure H2O is injected at 0.25
PV/day and an equal reservoir-volume rate is produced, for 1 PVI total.

The current overall composition is H2O/squalane = 0.60/0.40 mole fraction.
It is initialized by a P-T-z flash as one hydrocarbon-rich phase with exactly
zero water saturation.  Pure-water injection must therefore trigger a genuine
phase-appearance event; no finite initial water phase is prescribed.

The case answers two code-level questions before adding more pseudo-components:

1. Can the production H2O/squalane PR flash remain coupled to Darcy flow while
   the SCW front creates and removes phases?
2. Can the simulator report equilibrium phase viscosities and the producer's
   signed H2O/heavy-component mass rates without reconstructing them from a
   cell sample?

The squalane `Tc/Pc/omega/kij` match the existing Teratani-based flash input.
Its critical volume is the same PR-derived placeholder used by that tool.
Critical volume does not affect the unshifted-PR fugacity calculation, but it
does affect LBC viscosity.  Consequently phase equilibrium is tied to the
existing flash baseline, while absolute heavy-phase viscosity is screening
level until `Vc` or a viscosity multiplier is calibrated to data.

Run from the repository root (WSL/Linux profile):

```bash
make case CASE=scw_kerogen_squalane_1d -j
make run CASE=scw_kerogen_squalane_1d NP=1 \
  RESULT_DIR=results/base
python case/scw_kerogen_common/analyze_results.py \
  case/scw_kerogen_squalane_1d/results/base
```

Primary outputs are `reservoir_diagnostics.csv` (oil/gas/water viscosity
min/mean/max), `well_history.csv` (signed phase and component mass rates),
`solution_*.csv`, and `phase_state_*.csv`.  The analyzer creates
`producer_composition.csv`, `viscosity_history.csv`, `summary.json`, and a PNG
when Matplotlib is available.

The previously checked 1-PVI run below belongs to the superseded 0.75/0.25,
`kij=0.2395` baseline and is retained only as historical evidence.  It reached 4 days with maximum whole-run
component-balance relative error `6.32e-8`.  Its final volume-weighted viscosities are about
`9.36e-4 Pa s` (oil-rich slot) and `6.62e-5 Pa s` (water-rich slot); the
producer H2O mass fraction about 0.853.  It must not be combined with the
current physical-baseline parameters.

At 653.2 K the hydrocarbon-enriched dense fluid is not guaranteed to occupy the
simulator's conventional oil slot.  Interpret `phase_presence_mask`,
saturations, and all O/G/W viscosities together; do not relabel a gas-slot
SCW-rich dense phase as conventional gas or oil solely from the slot name.

The full staged design and recommended sensitivity matrix are documented in
[`../../docs/SCW_KEROGEN_FLOW_EXPERIMENT.md`](../../docs/SCW_KEROGEN_FLOW_EXPERIMENT.md).
