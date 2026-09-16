# three_eos_3d_compare

## 1. Purpose

This is the project's common **3-D structured-grid PR / Søreide–Whitson / CPA comparison case**. The three runs use the same reservoir, wells, initial pressure-temperature-overall-composition state, transport discretization, relative permeability, component set and non-aqueous binary-interaction baseline. The executable selects the thermodynamic backend with `-eos pr|sw|cpa`.

The case intentionally uses only the new `FullyCompositionalThreePhase` model. H2O, CO2 and all hydrocarbons can partition among the oil-rich, gas-rich and water-rich phases through the common fugacity-equilibrium flash. The following legacy optional models are disabled:

- old independent aqueous-CO2 dissolution variable;
- adsorption;
- Land trapped-gas model.

The case is a controlled EOS-comparison benchmark, not a history-matched field model.

## 2. Components and why this set was chosen

The common five-component fluid is

```text
H2O / CO2 / CH4 / C2H6 / nC4H10
```

n-butane is deliberately the heaviest component. The original Søreide–Whitson hydrocarbon-aqueous correlation was fitted over CH4 through nC4, so this case does not require the heavy-hydrocarbon extrapolation used by some reservoir pseudo-component examples.

At the common initial state

```text
T = 305 K
P = 60 bar
z = [0.30, 0.10, 0.15, 0.15, 0.30]
```

all three production flash backends resolve O+G+W. Current regression values are:

| EOS | So | Sg | Sw |
|---|---:|---:|---:|
| PR | 0.59379051 | 0.33291295 | 0.07329653 |
| SW | 0.59499883 | 0.33094096 | 0.07406021 |
| CPA | 0.56264318 | 0.37788412 | 0.05947270 |

These values are outputs of the P-T-z flash, not user-prescribed saturations.

## 3. Pure-component data

Critical properties, acentric factors and molar masses are taken from the public CoolProp 8.0 fluid-property pages. Critical volume is `1/rho_c,molar` from the same source.

| component | Tc [K] | Pc [MPa] | Vc [m3/mol] | omega | M [kg/mol] |
|---|---:|---:|---:|---:|---:|
| H2O | 647.096 | 22.0640 | 5.594803743e-5 | 0.3442920843 | 0.018015268 |
| CO2 | 304.1282 | 7.377298 | 9.411848339e-5 | 0.22394 | 0.0440098 |
| CH4 | 190.5640 | 4.599200 | 9.862771707e-5 | 0.01142 | 0.0160428 |
| C2H6 | 305.322 | 4.8722 | 1.458387816e-4 | 0.099 | 0.03006904 |
| nC4H10 | 425.125 | 3.796000 | 2.549219298e-4 | 0.200810094644 | 0.0581222 |

Public software source:

- CoolProp fluid documentation, `https://coolprop.org/fluid_properties/fluids/`.

The underlying CoolProp pages cite their pure-fluid reference EOS data sources separately.

## 4. PR parameters

PR uses the original Peng–Robinson cubic form implemented by the project:

- `OmegaA = 0.45724`;
- `OmegaB = 0.07780`;
- original PR kappa branch (`eosModelFlag = 1`);
- classical quadratic mixing rule.

Reference:

- Peng, D.-Y.; Robinson, D. B. *A New Two-Constant Equation of State*, Industrial & Engineering Chemistry Fundamentals 15 (1976) 59-64. DOI `10.1021/i160057a011`.

## 5. Common non-aqueous BIP baseline

To keep the three reservoir runs comparable, all three backends use the same non-aqueous `k_ij` matrix. It is not fitted to the synthetic geology.

```text
              H2O      CO2      CH4      C2       nC4
H2O           0       0.1896   0.4850   0.5000   0.5000
CO2          0.1896    0       0.1000   0.1300   0.1277
CH4          0.4850   0.1000    0        0       0.09281
C2           0.5000   0.1300    0        0        0
nC4          0.5000   0.1277   0.09281   0        0
```

The H2O row follows the non-aqueous water-interaction pattern already used by the project's SW benchmarks. The CO2/light-hydrocarbon entries follow the published SPE3-style interaction pattern used in the project's validated Panfili/SPE3 closure. Unspecified light-hydrocarbon pairs remain zero rather than being tuned for this case.

This is a **comparison convention**: keeping this matrix common makes the observed reservoir differences easier to attribute to PR vs SW water treatment vs CPA association rather than to independent per-EOS BIP fitting.

Parameter traceability for this common matrix:

- H2O-CO2 / H2O-CH4 / H2O-heavier-hydrocarbon non-aqueous values follow the Søreide-Whitson non-aqueous pattern used by the project's SW validation path;
- the CO2-C2 (`0.13`), CO2-C4-6 (`0.1277`) and C1-C4-6 (`0.09281`) values are reproduced in the published SPE3 interaction table; nC4 is represented by the C4-C6 entry;
- Raimondi, L. *Comment on “Three-Phase Equilibrium Computations for Hydrocarbon-Water Mixtures Using a Reduced Variables Method”*, Ind. Eng. Chem. Res. 59 (2020) 3287-3289, DOI `10.1021/acs.iecr.9b05598`, reproduces the SPE3 BIP table used here for traceability;
- zero entries are explicit benchmark choices, not fitted values.

## 6. Søreide–Whitson parameters

SW uses the same non-aqueous matrix above and adds the original Søreide–Whitson water treatment:

- water-specific PR alpha function;
- aqueous H2O-CO2 correlation;
- aqueous H2O-hydrocarbon correlations for CH4, C2H6 and nC4H10;
- fresh water, salinity = 0 mol/kg H2O.

Reference:

- Søreide, I.; Whitson, C. H. *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with pure water and NaCl brine*, Fluid Phase Equilibria 77 (1992) 217-240. DOI `10.1016/0378-3812(92)85105-H`.

The choice of nC4 as the heaviest hydrocarbon keeps every hydrocarbon inside the stated fit range of the implemented SW hydrocarbon-aqueous correlation.

## 7. CPA parameters

CPA uses the project's SRK-CPA implementation with explicit pure-component parameters. The table below is transcribed from Table 1 of Qvistgaard et al. (2023). The published table gives `b`, `Gamma` and `c1`; the code uses the equivalent project convention `a0 = b R Gamma`.

| component | b [cm3/mol] | Gamma [K] | c1 | association |
|---|---:|---:|---:|---|
| H2O | 14.52 | 1017.3 | 0.6736 | 4C, beta=0.0692, epsilon/R=2003.2 K |
| CO2 | 27.2 | 1551.22 | 0.7602 | none |
| CH4 | 29.1 | 959.02 | 0.4472 | none |
| C2H6 | 42.9 | 1544.54 | 0.5846 | none |
| nC4H10 | 72.081 | 2193.08 | 0.7077 | none |

References:

- Qvistgaard, D. et al. *Parameterization and Uncertainty Analysis of Binary Interaction Parameters for Triethylene Glycol and Ethane/Propane*, Fluid Phase Equilibria 570 (2023) 113796. DOI `10.1016/j.fluid.2023.113796`.
- Yakoumis, I. V. et al. *Prediction of Phase Equilibria in Binary Aqueous Systems Containing Alkanes, Cycloalkanes, and Alkenes with the Cubic-plus-Association Equation of State*, Ind. Eng. Chem. Res. 37 (1998) 4175-4182. DOI `10.1021/ie970947i`. This work shows why 4C water is appropriate for aqueous alkane mutual-solubility calculations.
- Kontogeorgis, G. M. et al. original CPA framework: Fluid Phase Equilibria 118 (1996) 27-59, DOI `10.1016/0378-3812(95)02843-9`.

Only H2O self-associates in this five-component case. No case-specific CPA parameter is fitted to the resulting reservoir output.

## 8. Geology

The model is deterministic and synthetic so every EOS sees exactly the same spatial problem.

```text
grid      = 36 x 24 x 8 = 6912 cells
domain    = 1200 x 720 x 96 m
cell size ~ 33.3 x 30 x 12 m
```

Eight layers use horizontal permeability between 5 and 220 mD and porosity between 0.12 and 0.225. Layer `k=3` is a low-permeability baffle. A localized window at `i=16..20` along the diagonal corridor raises the baffle permeability to approximately 80/65/6 mD in x/y/z, so cross-layer flow is focused through a finite opening rather than uniformly leaking everywhere. A five-cell-wide diagonal channel connects the lower southwest region to the upper northeast region; outside the baffle it doubles horizontal permeability and increases porosity by 0.015.

This range represents a moderately heterogeneous sandstone model. It is an engineering benchmark choice, not an attempt to reproduce one named field.

## 9. Wells

Two vertical wells force a genuinely three-dimensional displacement path:

- `CO2_INJ`: `(i,j)=(3,4)`, completed in `k=0..2`; pure CO2 gas-slot injection, 300000 surface/reference m3/day;
- `PROD`: `(i,j)=(32,19)`, completed in `k=4..7`; BHP = 50 bar.

The gas surface/reference density is 1.80 kg/m3, close to a standard-condition CO2 reference density; oil and water reference densities are 620 and 998 kg/m3. The rate-control equation in the current well model uses **surface phase volume rate**, so the injector target above is intentionally reported on that basis.

The initial reservoir pressure is 60 bar. The injector and producer lie on opposite sides of the low-permeability baffle and near the diagonal channel, so horizontal displacement, vertical pressure communication and phase redistribution all matter.

## 10. Time control

Default simulation:

```text
120 common output intervals x 0.5 day = 60 days
console/diagnostic output every 1 common interval = 0.5 day
internal adaptive stepping = ON
```

The 0.5-day schedule is the **common physical output grid**, not a requirement that every accepted backward-Euler step must be 0.5 day. PR, SW and CPA use the same adaptive policy. If one EOS cannot converge over a full 0.5-day attempt, the transaction is rolled back, the internal step is cut, and the stepper still lands exactly on the next common 0.5-day output time. Therefore reservoir and well histories remain directly comparable at identical physical times.

Different accepted/rejected internal-step counts are retained deliberately: they quantify nonlinear robustness and computational cost of each EOS and should be reported separately from the physical-output comparison.

From v59, phase appearance/disappearance still uses the single v58 hysteretic active set, but `S<=1e-4` is explicitly treated as a **thermodynamic probe band**, not a physical residual saturation: it only opens a restricted reduced-phase candidate and missing-phase stability remains the final guard. A phase removed from this trace band records suppression history and uses a `1e-4` reappearance deadband; naturally absent phases keep the stricter `1e-6` appearance criterion. The fully-compositional Newton limiter is also phase-aware, so a pathological composition correction in a disappearing phase cannot damp every other phase and all saturations. `phase_presence_mask` keeps the historical 1..7 encoding, while `phase_hysteresis_suppression_mask` exposes the numerical history without changing the physical output grid.

This case also enables a residual-stagnation guard for the expensive failed attempts seen in the v58 PR/SW/CPA logs. PETSc's normal convergence test remains authoritative; only a persistent plateau after at least 8 Newton iterations, with 5 iterations lacking `1e-4` relative improvement, is terminated early so the existing transactional rollback and dt cut can proceed. `SNES max_it=50` is deliberately retained for difficult solves whose residual continues to improve.

The current production case uses the single maintained Natural/active-set path documented in `docs/MODEL.md` and `docs/THERMODYNAMICS.md`; no historical experimental phase-boundary/globalization runtime switches are part of this case.

The production nonlinear path is intentionally small: a single hysteretic active set, the `1e-4` thermodynamic probe followed by restricted Flash and missing-phase stability, phase-aware Newton limiting, transactional rollback, and the residual-stagnation guard. PR/SW/CPA share the certified reduced-set fallback and strict stability-failure semantics documented in `docs/THERMODYNAMICS.md`.

Consistent Natural residual/variable scaling is retained and enabled by default because it changes coordinates, not the physical zero of the equations. It may be disabled for a clean scaling A/B:

```bash
EOS=cpa sbatch run.sh -natural_scaling 0
```

Scaling controls include `-natural_pressure_scale`, `-natural_composition_scale`, `-natural_saturation_scale`, `-natural_mass_residual_scale`, `-natural_fugacity_residual_scale`, `-natural_closure_residual_scale` and `-natural_rate_residual_floor`.

## 11. Build

From `case/`:

```bash
make three_eos_3d_compare -j
make prepare CASE=three_eos_3d_compare TASKS=8
```

The default launcher profile is 8 MPI ranks on `ft_module`.

## 12. Run one EOS

`run.sh.in` passes the `EOS` environment variable to `-eos`.

```bash
cd three_eos_3d_compare

EOS=pr  RESULT_DIR=./results/pr  sbatch run.sh
EOS=sw  RESULT_DIR=./results/sw  sbatch run.sh
EOS=cpa RESULT_DIR=./results/cpa sbatch run.sh
```

Equivalent direct PETSc option:

```bash
./bin/three_eos_3d_compare -eos sw -result_dir ./results/sw ...
```

If no `-eos` is given, PR is the default. If no result directory is given, the executable automatically uses `results/pr`, `results/sw` or `results/cpa` according to the selected backend.

## 13. Submit all three runs

After compiling and preparing `run.sh`:

```bash
cd three_eos_3d_compare
./submit_all_eos.sh
```

This submits three independent Slurm jobs with identical numerical settings and separate result directories.

## 14. Compare output

After all three jobs finish:

```bash
python3 compare_eos_results.py
python3 plot_eos_results.py --force
```

It creates:

```text
results/comparison/eos_reservoir_history.csv
results/comparison/eos_well_history.csv
results/comparison/eos_solver_summary.csv
results/comparison/eos_final_summary.csv
results/figures/{en,zh}/*.png
results/figures/{en,zh}/*.pdf
results/figures/figure_manifest.json
```

`plot_eos_results.py` creates eight separate, single-axis figures in both
English and Chinese. Legends remain English in both versions. The figures have
no titles or grid lines, retain a four-sided frame, and distinguish PR/SW/CPA
with color, line style and marker shape. PNG files are exported at 120 x 80 mm
and 300 dpi; PDF files retain vector text and line art.

The plotting script reads every common output time without smoothing,
interpolation or resampling. Producer rate is displayed as positive withdrawal
magnitude (`-q_total_surface`). The mass-balance figure shows the maximum
absolute relative error across the five components at each output time and
uses a symmetric-logarithmic axis so the exact-zero initial value is retained.
The generated manifest records source-file hashes, transformations, styles,
dimensions, software version and output hashes.

| Figure stem | Quantity | Scientific reading |
|---|---|---|
| `01_average_pressure` | volume-weighted average reservoir pressure | compares pressure support and depletion response |
| `02_average_oil_saturation` | average oil-rich-phase saturation | compares liquid-hydrocarbon redistribution |
| `03_average_gas_saturation` | average gas-rich-phase saturation | compares gas development during CO2 injection |
| `04_average_water_saturation` | average water-rich-phase saturation | isolates the largest aqueous-model difference |
| `05_producer_total_rate` | positive production-rate magnitude | compares the common 50-bar producer response |
| `06_injector_bhp` | CO2-injector bottom-hole pressure | compares pressure required to maintain the common rate target |
| `07_component_mass_balance` | maximum absolute component-relative error | verifies conservation without hiding the exact-zero initial state |
| `08_time_step_attempts` | accepted plus rejected internal attempts | compares nonlinear robustness separately from physical output |

The underlying accessible data are the CSV files in `results/comparison/` and
the component mass-balance files in each EOS result directory. Wall-clock time
is deliberately not plotted because scheduler placement and concurrent CPU
sharing can change it independently of the equations of state.

The most useful comparison quantities are:

- average pressure;
- average So/Sg/Sw;
- phase density and viscosity;
- O/G/W phase-presence cell counts;
- producer BHP and phase/total surface rates;
- component mass balance, especially injected/produced CO2;
- accepted/rejected internal steps, minimum accepted dt, attempted SNES/KSP work and simulation-loop wall time;
- detailed SNES/KSP histories in each EOS result directory.

`eos_solver_summary.csv` keeps numerical robustness/cost separate from the physical reservoir comparison. Do not interpret different wall-clock time alone as an EOS-accuracy ranking. Compare physical outputs and nonlinear/linear solver cost separately.

## 15. Validation guard

`test/src/unit/three_eos_3d_case_test.cpp` uses the same production `fluid_factory.hpp` and production three-phase flash to verify:

- all legacy optional physics are disabled;
- grid is genuinely three-dimensional;
- common P-T-z is normalized;
- PR, SW and CPA all converge to O+G+W at the common initial state;
- initial saturation anchors remain unchanged unless a deliberate thermodynamic change is made.

## Failed mass-balance / near-well flux diagnostics

From v52 onward the benchmark keeps adaptive `[TIME][TARGET]`, `[TIME][RETRY]`, `[TIME][ACCEPT]` and `[TIME][REJECT]` events enabled and prints every common output target by default. A nonlinear failure is diagnosed **before** rollback. The old phase-transition diagnostic has been removed because target-cluster v51 runs showed no phase-presence changes at the residual plateau.

Each failure now locates the MPI-global largest **component mass-conservation** residual and prints:

- diagnostic cell/current-input id, target component, volume and pressure;
- current phase saturation, density, mobility and phase mass fractions used by transport;
- for every component, `accumulation + face_flux - well_source`, reconstructed residual, PETSc residual and reconstruction error;
- for every neighboring TPFA face and every phase: potential difference, upwind side, mobility, density, Darcy rate, phase mass rate and per-component mass flux;
- for every perforation on the diagnostic cell: well/control/target/actual/control residual/BHP/WI, phase pressure drop, mobility, density, surface/reservoir rate, phase mass rate, injection fractions and component source.

Machine-readable diagnostics are written to:

```text
results/<eos>/failed_mass_balance_diagnostics.csv
results/<eos>/failed_face_flux_diagnostics.csv
results/<eos>/failed_well_source_diagnostics.csv
```

This instrumentation is diagnostic only. It does not change EOS, flash, phase-state criteria, TPFA/upwind physics, well equations, residual/Jacobian or adaptive time-step policy.

For CPA, v57 also writes cumulative thermodynamic profiling to:

```text
results/cpa/thermodynamic_profile.csv
```

The profiler reports MPI-summed call counts and max/sum-rank wall times for CPA phase evaluations, mixing/cache hits, association work and density-root work. Timings are nested rather than exclusive. The water-only 4C case uses an analytic site-fraction fast path and an isothermal temperature-coefficient cache; general cross-association still uses the original iterative path.

## 16. MRST geology preview

`scripts/plot_geology_mrst.m` builds a deliberately simple proposed `20 x 20 x 5`
geology with the local MRST installation and exports three separate views in
English and Chinese. The model contains background sandstone, one straight
high-permeability channel, one middle shale layer and one `4 x 4` flow window:

- a vertically exaggerated three-dimensional overview;
- a top view of the reservoir channel;
- a top view of the shale baffle and its flow window.

Run it from MATLAB with:

```matlab
run('case/three_eos_3d_compare/scripts/plot_geology_mrst.m')
```

The figures, cell-wise CSV data and transformation manifest are written below
`results/geology/`. The geometry is deterministic; the script applies no
interpolation, smoothing or stochastic realization. This is a design preview
only: it does not change the production grid or reservoir solver configuration.
The presentation figures retain only the geology, wells and in-frame legend;
coordinate axes, ticks, labels and internal cell-edge grids are hidden. The 3-D
wellbores extend from above the model to the full completion interval. A dark
outer frame retains the full geological-grid boundary in all three views.
