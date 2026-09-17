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

- **Locked now:** the experimentally observed four boiling-range boundaries are the preferred lumping topology for the case.
- **Directly measured for paired pure kerogen:** the 380 °C recovered-oil mass fractions `0.81 / 23.73 / 34.11 / 41.35 wt%`.
- **Not yet locked for intact raw shale:** the four lump mass fractions for the RSC raw-shale product, because a same-physical-sample simulated-distillation table has not yet been recovered.

Until same-sample raw-shale distillation is found, the pure-kerogen fractions are an explicit paired experimental prior, not silently relabelled primary data.

## No representative molecule is assigned at the lumping stage

Lumping defines **what material range belongs to each pseudo-component**. It does not define thermodynamic properties.

The following remain unset until supported by data/correlations appropriate to each measured fraction:

- average molecular weight;
- density / specific gravity;
- mean or characterization boiling temperature;
- critical temperature and pressure;
- acentric factor;
- PR/SW binary interaction coefficients;
- CPA association scheme and parameters;
- viscosity-correlation parameters.

In particular, `OIL_HEAVY` must not be equated to squalane merely because squalane is computationally convenient. The experimental SARA evidence shows substantial resin/asphaltene content, so a nonpolar single-normal-alkane surrogate can only be introduced as a documented approximation and must be validated independently.

## Property-characterization sequence

For each experimental boiling-range lump, parameterization must proceed in this order:

1. recover or estimate `MW` and `SG/density` from experimental distillation/GC/product characterization;
2. define a representative boiling temperature or full boiling-range characterization;
3. derive/check `Tc`, `Pc`, `omega` using a documented petroleum-fraction correlation and uncertainty range;
4. calibrate H2O-lump and lump-lump interaction parameters against high-T/P VLE/LLE/PVT data;
5. validate phase behavior at `380 °C / 25 MPa` before any flow run;
6. separately validate viscosity/transport closure.

No step may be replaced by selecting a convenient normal alkane and declaring it to be the experimental fraction.

## Acceptance criterion

The lumping stage is considered complete when:

- all case documentation uses `OIL_GASOLINE / OIL_DIESEL / OIL_MIDDLE / OIL_HEAVY` or the corresponding boiling ranges;
- no production-case composition is defined by `nC4 / nC10 / squalane` or arbitrary carbon-number cuts;
- experimental mass fractions and their sample/evidence status are stored explicitly;
- missing EOS/transport properties remain explicit gaps rather than filled with screening surrogates.
