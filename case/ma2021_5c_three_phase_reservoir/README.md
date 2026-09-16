# Ma et al. (2021) five-component O/G/W mutual-solubility benchmark

This case replaces the previous `yang2019_1d_co2_flood` case because the Yang
model keeps water as an independent immobile phase and only transfers CO2 into
water.  The present case is selected specifically to exercise the **fully
compositional O/G/W model in which H2O and hydrocarbons belong to the same
three-phase equilibrium calculation**.

## Literature source

X. Ma et al., *Three-Phase Equilibrium Calculations of
Water/Hydrocarbon/Nonhydrocarbon Systems Based on the Equation of State (EOS)
in Thermal Processes*, **ACS Omega** 6(50) (2021) 34406-34415.
DOI: `10.1021/acsomega.1c04522`.

The input comes from Tables 7-8 and the reference split from Table 9:

- `P = 13.79 bar`, `T = 366.5 K`;
- overall mole fraction `[H2O,C1,C6,C10,C15] = [0.10,0.10,0.20,0.40,0.20]`;
- Table-7 critical properties and acentric factors (the PMC text rendering
  shows `C15 Pc=08.49 bar`; this is treated as a transcription defect and
  `18.49 bar` is used, consistent with the same C15 property in Table 4);
- Table-8 **non-aqueous** PR binary-interaction coefficients;
- Table-9 reference phase compositions and phase mole fractions.

The published result is genuinely three phase.  For example, H2O is about
0.622 mol% in the oil-rich phase and 5.9951 mol% in the gas-rich phase, while
C1 is finite (about 0.0004 mol%) in the water-rich phase.  Thus this is not a
free-water or Henry-only model.

## Important model boundary

Ma et al. use a modified Peng-Robinson/Soreide-Whitson treatment.  The aqueous
phase and the non-aqueous phases use different H2O-related BIP sets, and the
water alpha function is modified.  **v24 now provides the original selectable
Soreide-Whitson backend**, but Ma et al. use their own modified aqueous
parameterization rather than the untouched original SW correlations.  This case
therefore remains explicitly on ordinary PR until the exact Ma parameter set is
implemented and validated; silently switching it to original SW would not be an
exact literature reproduction.

Therefore this case is intentionally split into two validation levels:

1. **Exact literature input/reference data** are stored in `case_config.hpp` and
   checked by `ma2021_case_test`.
2. **MPMC ordinary-PR result** must still produce stable O+G+W equilibrium and
   true mutual partitioning.  It is compared with Table 9, but exact equality of
   phase fractions is *not* asserted until a phase-dependent aqueous PR/SW
   backend is implemented.

This distinction prevents the project from claiming a quantitative
reproduction that the current thermodynamic model cannot yet support.

## Reservoir-flow wrapper

The paper is a phase-equilibrium benchmark and does not define a porous-medium
flow geometry or wells.  The following are therefore **MPMC test parameters,
not literature parameters**:

- structured grid `40 x 1 x 1`, domain `200 x 10 x 10 m`;
- porosity `0.20`, permeability `100/100/10 mD`;
- pure-C1 gas injector at `1 std m3/day` on the left;
- producer at `13.50 bar` on the right;
- isothermal `366.5 K`, 20 output intervals of `0.1 day`.

The purpose of this wrapper is to run the literature three-phase fluid through
PETSc residual/Jacobian assembly, Darcy transport, wells, adaptive stepping and
`component_mass_balance.csv` without pretending these flow parameters came from
Ma et al.

## Build and run on the PETSc cluster

```bash
cd case
make ma2021_5c_three_phase_reservoir -j
make reset-run CASE=ma2021_5c_three_phase_reservoir
cd ma2021_5c_three_phase_reservoir
./run.sh
```

For a very short first smoke run:

```bash
./bin/ma2021_5c_three_phase_reservoir -numSteps 2 -dt 0.05 -out_step 1
```

Inspect especially:

- the initial `Phase-state cells` line (should start as `OGW` for the current
  ordinary-PR comparison);
- Newton/KSP convergence;
- `results/component_mass_balance.csv` for H2O, C1, C6, C10 and C15.
