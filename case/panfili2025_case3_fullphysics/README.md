# Panfili et al. (2025) use case #3 — full-physics live-water CO2 storage

This case is the first MPMC_SCW production benchmark derived from a published
**multicomponent oil/gas/water reservoir simulation using the Soreide-Whitson
EOS with all primary components allowed in all three reservoir phases**.

The selected target is Panfili et al. (2025), use case #3, with the Section 7.3
variant requested for comparison:

- depleted gas-condensate reservoir;
- 30% connate-water saturation;
- pure CO2 injection;
- live-water/full-physics equilibrium: aqueous dissolution **and** water
  vaporization enabled;
- isothermal reservoir at 200 F;
- full O/G/W equilibrium with Soreide-Whitson thermodynamics.

## Literature source

P. Panfili, L. Patacchini, A. Ferrari, K. Esler and A. Cominelli,
“Three-phase equilibrium in a GPU-based compositional reservoir simulator,”
*Computational Geosciences* 29 (2025), article 31.
DOI: `10.1007/s10596-025-10369-3`.

The paper states that its live-water formulation allows every primary
component, including H2O, to occur in vapor, liquid and aqueous reservoir
phases and uses a full three-phase flash with the Soreide-Whitson EOS.

## What is taken directly from the paper

### Reservoir / rock summary

- grid dimensions: `47 x 99 x 15`;
- average cell dimensions: `328 x 328 x 8 ft`;
- average horizontal permeability: `815 mD`;
- average porosity: `15%`;
- initial pressure: `3000 psi` at reference/GWC depth `2910 ft`;
- reservoir temperature: `200 F`;
- Section-7.3 connate-water case: `Swc = 30%`.

**Paper inconsistency retained in the documentation:** `47 x 99 x 15` equals
69,795 cells, whereas the same paragraph says “approximately 400k active
cells”.  The runnable structured surrogate uses the explicit dimensions rather
than silently changing them.

### Panfili Table 3 non-water fluid

The nine non-water mole fractions and properties are stored verbatim (with
unit conversion to SI) in `case_config.hpp`:

| component | mole fraction | MW | Tc (R) | Pc (psia) | omega |
|---|---:|---:|---:|---:|---:|
| CO2 | 0.01210 | 44.01 | 548.46 | 1071.33 | 0.22500 |
| N2 | 0.01940 | 28.013 | 227.16 | 492.31 | 0.04000 |
| C1 | 0.65990 | 16.043 | 343.08 | 667.78 | 0.01300 |
| C2 | 0.08690 | 30.07 | 549.77 | 708.34 | 0.09860 |
| C3 | 0.05910 | 44.097 | 665.64 | 618.70 | 0.15240 |
| C4-6 | 0.12970 | 66.869 | 806.54 | 514.93 | 0.21575 |
| C7+1 | 0.02745 | 107.779 | 838.11 | 410.75 | 0.31230 |
| C7+2 | 0.00515 | 198.562 | 1058.04 | 247.56 | 0.55670 |
| C7+3 | 0.00030 | 335.198 | 1291.89 | 160.42 | 0.91692 |

H2O is the tenth component and uses the project's validated SW water model.

### Well schedule

Paper-defined schedule retained in the wrapper:

1. years 0-32: two crest producers, each `50,000 Mscf/day` surface gas rate,
   minimum BHP `435 psi`;
2. years 32-34: idle period;
3. from year 34: the two wells become injectors, each `90,000 Mscf/day`,
   maximum BHP `2850 psi`;
4. Section 7.3 changes the injection stream to **pure CO2** and enables both
   dissolution and water vaporization.

The implementation keeps **exactly two physical `Well` objects**, matching the
paper statement that the same two crest wells are converted from producers to
injectors.  Their operating mode is updated at each attempted time-step end:
`Depletion -> Idle -> Injection -> Closed`.  This is important because Natural
stores one implicit BHP unknown/equation at each well representative cell; the
old v25 staging workaround (separate producer and injector objects on the same
cell) was invalid and was rejected before SNES setup.

For the one-year backward-Euler comparison grid, the schedule convention is:
- the `31 -> 32 y` step is still a production step;
- `32 -> 34 y` is idle;
- the first injection step is `34 -> 35 y`;
- `49 -> 50 y` is the final figure-derived injection step.

Figure 25 shows cumulative injection becoming flat at approximately year 50,
but the text does not state the exact injection-stop time.  The default wrapper
therefore uses year 50 as an explicitly **figure-derived comparison value**;
it is not presented as an exact paper schedule datum.

The default simulation horizon is 100 years to match the time scale of Figs.
25-26.  The paper also discusses >1000-year plume monitoring; users can extend
`Time::numberOfSteps` after the short comparison run is validated.

## Thermodynamically consistent 30% connate-water initialization

The paper gives **water saturation**, not overall H2O mole fraction.  Setting
`z_H2O=0.30` would therefore be wrong because saturation and phase mole fraction
are not equal.

The case preserves the Table-3 ratios of all nine non-water components and
solves the production SW P-T-z flash at 3000 psi / 200 F for the overall H2O
fraction that gives `Sw=0.30`.  The resulting overall composition is stored in
`InitialState::overallComposition`.

The numerical reference is also stored in `reference_initial_flash.csv`; the pressure sweep is in `reference_pressure_sweep.csv`.

Current production flash result:

- phase state: `G+W`;
- `Sg = 0.700000`, `Sw = 0.300000`, `So = 0`;
- water is finite in the gas phase;
- CO2 and hydrocarbon components are finite in the water-rich phase.

A pressure sweep with the same feed gives `G+W` near 2000 psi and `O+G+W` by
about 1800 psi.  This protects the intended condensate appearance during
primary depletion.

## Data gaps: literature values vs. MPMC closures

This repository deliberately does **not** claim a cell-for-cell reproduction of
the Hamilton field model.  The paper cites an open Hamilton CCS dataset, but the
archive itself is not bundled here.  Several inputs required by the current
runtime are not printed numerically in the paper.

| item | source/status in this runnable case |
|---|---|
| 47x99x15 dimensions, average cell dimensions | paper |
| average kh, porosity | paper |
| full heterogeneous porosity/permeability fields | Hamilton archive; **not bundled** |
| vertical permeability | not printed; wrapper uses `kz=815 mD` as explicit isotropic closure |
| exact producer/injector completion cells | not printed; two symmetric upper-layer cells are used |
| Table-3 non-water PVT properties | paper |
| exact Panfili hydrocarbon-hydrocarbon BIPs | not printed; published SPE3 interaction pattern is used as an explicit closure |
| critical volumes for pseudo-components | not printed; `Zc=0.27` estimate used only by current LBC viscosity path |
| case-#3 salinity | not stated; wrapper explicitly assumes fresh water, `0 mol/kg`, pending source-deck recovery |
| relative permeability numeric tables | Fig.17 only; simple endpoint-aware Corey curves are isolated in `case_fluid.hpp` |
| aqueous density/viscosity model | Panfili uses brine correlations + Ezrokhi; current full-3p MPMC still uses its common EOS/LBC property path |
| surface separator train | Panfili has dehydration + separator stages; current MPMC well conversion is simpler |

These closures are intentional and centralized so they can be replaced without
changing the reservoir equations or SW flash when the exact Hamilton inputs are
available.

## Why the flash core changed in v25

The 10-component SW fluid exposed a numerical issue in the existing two-phase
Rachford-Rice helper.  Heavy components in a water-rich phase can have
`K << 1e-14`; evaluating `1 + beta*(K-1)` at `beta=1` suffers cancellation and
can become exactly zero even though the physical denominator is `K > 0`.

v25 evaluates the same expression in its stable form

```text
D = (1-beta) + beta*K
```

and only rejects genuinely non-positive/underflowed denominators.  SW water
pairs also receive a phase-aware fugacity-ratio K seed.  The ordinary-PR seed
path is unchanged.  `panfili2025_case3_test` protects this regression.

## Primary comparison outputs

For the default 0-100 year run, compare with paper Figs. 25-26 using:

- `results/component_mass_balance.csv` — CO2 inventory/injection/production and
  global mass-balance error;
- `results/reservoir_diagnostics.csv` — average reservoir pressure;
- `results/well_history.csv` — gas/water rates and BHP;
- `make_paper_comparison_csv.py` — combines those three standard files into
  `results/panfili_fig25_26_comparison.csv` with time in years, cumulative CO2
  in kg/Mton, average pressure in bar/psi and positive field-water-production
  magnitude in m3/s and m3/day;
- phase-state/solution snapshots — appearance of condensate during depletion and
  CO2 distributions after injection;
- `results/simulation_summary.csv` — SNES/KSP/runtime statistics.

For a future Hamilton-deck reproduction, Figs. 19-22 also provide useful
comparison targets: live/inert storage ratio versus connate Sw, CO2 mole
fractions in vapor/liquid/aqueous phases, and field oil-gas ratio.

## Build and run on the PETSc cluster

```bash
cd case
make panfili2025_case3_fullphysics -j
make reset-run CASE=panfili2025_case3_fullphysics TASKS=16
cd panfili2025_case3_fullphysics
./run.sh
```

A full 69,795-cell/10-component FIM run is not a small smoke case.  Before a
long run, use a short schedule override to verify the PETSc path and memory
requirements on the cluster.  Do not interpret a shortened schedule as a
literature comparison.
