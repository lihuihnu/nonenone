# SCW kerogen experiment-driven five-component 2-D screening

This case is the first 2-D flow screening for the current experiment-driven
five-component topology:

1. H2O
2. OIL_GASOLINE (IBP-180 C)
3. OIL_DIESEL (180-350 C)
4. OIL_MIDDLE (350-500 C)
5. OIL_HEAVY (>500 C)

It deliberately does **not** reuse the legacy nC4/nC10/squalane LMH identity.

## Scope

`CONDITIONAL_5C_FLOW_SCREENING_NOT_FORMAL_VALIDATION`

The purpose is to prove that the accepted five-component PVT kernel can be
connected to native 2-D flow, exact conserved-component producer output,
component-wise recovery, Light/Heavy enrichment and strict mass closure.

The run is not the formal 360/380 C paired laboratory claim. The formal
experiment remains blocked by apparatus/specimen measurements and independent
interaction/transport calibration gates.

## Frozen screening baseline

- T = 653.15 K (380 C)
- pressure center = 28 MPa
- grid = 60x20x1 over 0.30x0.10x0.010 m
- phi = 0.25
- k = 1e-12 m2
- initial z_H2O = 0.20
- oil mole ratio = characterized Gasoline/Diesel/Middle/Heavy ratio
- pure-water reservoir-rate injection = 2.5e-8 m3/s
- injector maximum BHP = 30 MPa
- producer BHP = 28 MPa
- effective PV used for actual-PVI accounting = 7.5e-5 m3
- PR76 five-component screening BIPs are frozen at 653.15 K
- HC-HC BIPs = 0

## Required outputs

The generic first-class producer writer must emit:

- RF_OIL_GASOLINE
- RF_OIL_DIESEL
- RF_OIL_MIDDLE
- RF_OIL_HEAVY
- RF_Light_group, with Light=Gasoline+Diesel
- RF_Heavy_group
- RF_total_hydrocarbon
- E_L_over_H_instant_mass
- E_L_over_H_cumulative_mass
- exact instantaneous/cumulative producer mass fractions
- phase-by-component producer ledgers

The trajectory must independently close H2O and all four hydrocarbon lumps.

## Interpretation

A successful trajectory establishes only numerical five-component flow
connectivity and conditional selectivity behavior under the current provisional
property/BIP package. It does not validate real Chang-7 selectivity magnitudes.
