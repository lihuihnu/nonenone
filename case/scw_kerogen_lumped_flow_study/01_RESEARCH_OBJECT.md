# 01 — Research object lock

## Decision

The primary geological precursor for this case collection is fixed to the **low-maturity, organic-rich lacustrine Chang 7 shale collected from outcrops in Tongchuan, southern Ordos Basin, China**, reported by Zhao et al. in *Sustainable Energy & Fuels* (2023), DOI: `10.1039/D2SE01361D`.

The target fluid for later flow simulations is **not generic kerogen**. It is the **pre-generated oil produced from this specific Chang 7 shale by supercritical-water conversion at 380 °C and 25 MPa**. This distinction must be preserved in all later PVT, lumping, viscosity and flow work.

## Source-rock identity

Reference: Qiuyang Zhao, Yu Dong, Lichen Zheng, Tian Xie, Baercheng Bawaa, Hui Jin and Liejin Guo, “Sub- and supercritical water conversion of organic-rich shale with low-maturity for oil and gas generation: using Chang 7 shale as an example,” *Sustainable Energy & Fuels*, 2023, 7, 155–163, DOI `10.1039/D2SE01361D`.

The reported sample is:

- Stratigraphy: Triassic Yanchang Formation, Chang 7 Member.
- Provenance: outcrop sample from Tongchuan, southern Ordos Basin, China.
- Depositional character: lacustrine shale.
- Kerogen type: Type II.
- Thermal maturity: low maturity, `Ro = 0.36–0.38 %`.
- TOC: `15.11 wt%`.
- Rock-Eval `S1 = 4.31 mg/g rock`.
- Rock-Eval `S2 = 76.35 mg/g rock`.
- Rock-Eval `Tmax = 443 °C`.
- Hydrogen index: `HI = 505.29 mg HC/g TOC`.
- Reported bulk elemental analysis: `C = 15.84 wt%`, `H = 2.33 wt%`, `N = 2.71 wt%`, `S = 3.80 wt%`.

The same paper identifies dolomite, feldspar, quartz, clay and pyrite as relevant mineral classes and shows that minerals materially affect oil yield under SCW conditions. Exact mineral fractions for the chosen sample remain to be extracted from the primary/supplementary data before any mineral-specific model is parameterized.

## Locked SCW conversion and numerical target condition

The direct experimental anchor and the intended primary numerical condition are now the same:

- Temperature: `380 °C` (`653.15 K`).
- Pressure: `25 MPa`.
- Duration of the reference conversion experiment: `4 h`.
- Water:shale mass ratio in the main temperature-series experiments: `1:1`.
- Reported oil yield at 380 °C: `352.1 mg/(g TOC)`.
- Oil characterization available in the paper: SARA fractions (saturates, aromatics, resins, asphaltenes).
- Gas characterization: gas chromatography; CH4, H2 and CO2 are explicitly discussed.

The paper notes that oil recovery by solvent evaporation can lose hydrocarbons with boiling points below about 46 °C; therefore the published liquid-oil composition must not be treated as a complete C1+ product distribution without correcting for this measurement boundary.

## Numerical-condition policy

The primary flow case shall use **380 °C / 25 MPa**. This removes the previous 25→28 MPa pressure extension and allows the thermodynamic and transport model to be compared directly against the closest available SCW conversion experiment before flow-specific assumptions are introduced.

The planned first temperature control is:

1. **Primary/reference flow state:** 380 °C, 25 MPa.
2. **Subcritical control:** 360 °C, 25 MPa, only after the fluid characterization is fixed.

Any later pressure sweep (for example 22–30 MPa) must be treated as a separate sensitivity study rather than part of the locked baseline.

## What is and is not locked by this decision

Locked:

- geological source: Tongchuan Chang 7 lacustrine shale;
- organic-matter class: low-maturity Type II kerogen;
- precursor geochemical baseline listed above;
- the 380 °C / 25 MPa SCW-generated oil from this sample as the fluid family to characterize;
- 380 °C / 25 MPa as the primary numerical thermodynamic condition;
- Zhao et al. (2023) as the primary source for the first-stage product characterization.

Not yet locked:

- Light/Middle/Heavy carbon-number cut points;
- gas lump definition;
- pseudo-component molecular weights, boiling points, specific gravities, critical properties or acentric factors;
- PR binary interaction coefficients;
- CPA pure-component/cross-association parameters;
- hydrocarbon-phase viscosity correlation;
- final initial overall composition for the flow case;
- permeability, porosity, relative-permeability and well-control values.

Those items belong to later research steps and must not be filled with screening values and then relabeled as experimental properties.

## Secondary literature — context only

A separate Chang 7 semi-open pyrolysis study reports a different sample with approximately `TOC = 26.8 wt%`, `Ro = 0.53 %`, Type II kerogen and `HI ≈ 476 mg/g TOC`, and shows strong generation of C6–C14, C15+ and wet gas around 360–380 °C. This literature is useful for defining candidate lump boundaries, but it is **not the same physical sample** as the Zhao et al. SCW sample and must not be mixed into the primary sample table without an explicit cross-sample label.

## Acceptance criterion for Step 1

Step 1 is complete only if all subsequent files and figures refer to the research object by the specific identity above instead of the generic phrases “kerogen” or “Chang 7 kerogen,” use 380 °C / 25 MPa as the locked primary condition, and label data imported from other Chang 7 samples with their own sample identifiers and evidence level.
