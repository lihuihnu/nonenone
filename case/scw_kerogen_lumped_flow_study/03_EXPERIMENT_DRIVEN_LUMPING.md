# 03 — Experiment-driven lumping

## Decision

The hydrocarbon pseudo-component structure for this case is defined **from the measured boiling-range separation of the generated oil**, not from preselected surrogate molecules or arbitrary carbon-number cuts.

The previous `nC4 / nC10 / squalane` and `C6-C14 / C15-C20 / C21+` ideas are no longer the default lumping definition for this research case. They may still be used in isolated code-regression or sensitivity tests, but they do not define the physical fluid.

## Experimental basis

Zhao et al., *Industrial & Engineering Chemistry Research* 2023, DOI `10.1021/acs.iecr.3c02759`, reports simulated-distillation bins for the generated oil. At `380 °C / 25 MPa`, Figure 6 directly labels:

| lump id | experimental name | boiling range | measured wt% of recovered oil |
|---|---|---|---:|
| OIL_GASOLINE | Gasoline | IBP-180 °C | 0.81 |
| OIL_DIESEL | Diesel | 180-350 °C | 23.73 |
| OIL_MIDDLE | Middle distillate | 350-500 °C | 34.11 |
| OIL_HEAVY | Heavy residue | >500 °C | 41.35 |

The four fractions sum to `100.00 wt%` on the **recovered generated-oil basis**.

These are not carbon-number bins and must not be silently renamed as `C6-C14`, `C15-C20`, `C21+`, or any individual hydrocarbon.

## Current fluid topology

For the first nonreactive SCW-flow model, the intended hydrocarbon representation is therefore:

1. `H2O`
2. `OIL_GASOLINE` — IBP-180 °C
3. `OIL_DIESEL` — 180-350 °C
4. `OIL_MIDDLE` — 350-500 °C
5. `OIL_HEAVY` — >500 °C

Gas products are **not merged into the oil lumps by default**. `H2`, `CO2`, `CH4`, `C2`, `C3`, `C4`, `C5`, ... may only be added after a common gas+liquid mass/mole basis has been reconstructed without double counting the low-boiling material lost during liquid recovery.

## Evidence status and sample scope

The 380 °C boiling-range fractions above come from the acid-pickled Type-II kerogen experiment and are therefore `SECONDARY_PAIRED` relative to the locked intact/raw-shale product dataset.

This distinction changes what can be locked now:

- **Locked:** the experimentally observed four boiling-range boundaries are the preferred lumping topology for the case.
- **Directly measured for paired pure kerogen:** the 380 °C recovered-oil mass fractions `0.81 / 23.73 / 34.11 / 41.35 wt%`.
- **Not yet locked for intact raw shale:** the four lump mass fractions for the RSC raw-shale product, because a same-physical-sample simulated-distillation table has not yet been recovered.

Until same-sample raw-shale distillation is found, the pure-kerogen fractions are an explicit paired experimental prior, not silently relabelled primary data.

## Characterization is separate from surrogate identity

Lumping defines **what material range belongs to each pseudo-component**. Property characterization then assigns an EOS starting property set to that measured range.

The first characterization pass is now documented in `04_PSEUDOCOMPONENT_CHARACTERIZATION.md` and `fluid_characterization/pseudo_component_characterization_380c.csv`. It uses:

- Figure 6 for lump mass fractions;
- Figure 5 for relative carbon-number structure inside each lump;
- petroleum SCN MW/Tb/SG characterization;
- Twu critical-property correlation;
- Edmister acentric-factor correlation.

This does not create a representative pure molecule. `representative_molecule` remains blank in the lump definition table.

In particular, `OIL_HEAVY` must not be equated to squalane merely because squalane is computationally convenient. Squalane is retained separately in `fluid_characterization/squalane_benchmark.csv` as a nonpolar heavy-saturated benchmark only.

## Current characterization status

The first-pass recovered-oil-basis properties are:

| lump | MW, g/mol | Tb, °C | SG | Tc, K | Pc, MPa | omega | Vc, cm3/mol |
|---|---:|---:|---:|---:|---:|---:|---:|
| `OIL_GASOLINE` | 82.392 | 68.991 | 0.691443 | 515.231 | 3.192179 | 0.269328 | 358.490 |
| `OIL_DIESEL` | 215.583 | 278.493 | 0.836212 | 736.423 | 1.751155 | 0.583473 | 824.073 |
| `OIL_MIDDLE` | 387.349 | 430.304 | 0.896458 | 870.437 | 1.138075 | 0.896554 | 1318.853 |
| `OIL_HEAVY` | 660.132 | 563.696 | 0.942605 | 982.864 | 0.808714 | 1.215676 | 1793.890 |

These are **provisional petroleum-characterization values**, not direct critical-property measurements. The Heavy row is explicitly `PROVISIONAL_HEAVY_TAIL_LOWER_BOUND_LIKE` because the GC-visible C38-C74 tail can under-represent resin/asphaltene material in the >500 °C fraction.

## Property-characterization sequence from here

1. retain the experimental four-bin topology and current provenance;
2. use the current table for initial PR/SW/CPA screening only;
3. validate or update MW/SG/Tb against same-fluid cut characterization when available;
4. calibrate H2O-lump and lump-lump interaction parameters against high-T/P VLE/LLE/PVT data;
5. treat Heavy association/PNA/SARA as a separate calibration problem rather than substituting squalane;
6. validate phase behavior at `380 °C / 25 MPa` before any flow run;
7. separately validate density/volume-shift and viscosity/transport closure.

No step may be replaced by selecting a convenient normal alkane and declaring it to be the experimental fraction.

## Acceptance criterion

The lumping/initial-characterization stage is considered complete when:

- all case documentation uses `OIL_GASOLINE / OIL_DIESEL / OIL_MIDDLE / OIL_HEAVY` or the corresponding boiling ranges;
- no production-case composition is defined by `nC4 / nC10 / squalane` or arbitrary carbon-number cuts;
- experimental mass fractions and their sample/evidence status are stored explicitly;
- the initial `MW/Tb/SG/Tc/Pc/omega/Vc` values are traceable to Figure 5/Figure 6 and documented characterization methods;
- squalane remains benchmark-only;
- unresolved high-T/P interaction and Heavy association data remain explicit gaps rather than being guessed.
