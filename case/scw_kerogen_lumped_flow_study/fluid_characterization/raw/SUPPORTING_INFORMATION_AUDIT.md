# Supporting-information audit

## Scope

This audit follows the locked research object in `01_RESEARCH_OBJECT.md`: the Tongchuan Chang 7 low-maturity raw shale studied by Zhao et al. at 380 °C, 25 MPa and 4 h. The purpose is to recover primary supplementary data before any Light/Middle/Heavy lumping.

## 1. RSC ESI — recovered

Primary article: Zhao et al., *Sustainable Energy & Fuels* (2023), DOI `10.1039/D2SE01361D`.

The electronic supplementary information contains Tables S1-S8. Tables S1-S3 resolve the main 380 °C raw-shale data gaps:

- Table S1: exact run conditions and charges for experiments 5 and 6 at 380 °C;
- Table S2: duplicate oil yields, SARA composition, residual shale data and mass balance;
- Table S3: duplicate total-gas and individual gas-component volumetric yields.

The exact duplicate values are preserved in `rsc_esi_380c_replicates.csv`. Means are allowed only as explicitly labelled derived summaries; raw replicates remain the authoritative observations.

### Important gas-basis consistency issue

The article text reports an H2 proportion of 26.9% at 380 °C. The ESI Table S3 reports component yields in mL/g shale. Directly dividing the tabulated H2 yield by the tabulated total-gas yield does not reproduce 26.9% for either 380 °C replicate. Therefore:

- the main-text 26.9% value is retained as a directly reported composition statement;
- the ESI mL/g-shale component values are retained as directly reported yield observations;
- no CH4/CO2/C2/C3+ percentages are invented by normalization until the reporting basis is reconciled;
- this discrepancy is a data-basis audit item, not a license to alter either source value.

## 2. ACS 2023 paired pure-kerogen study — main article recovered; detailed SI still outstanding

Paired study: Zhao et al., *Industrial & Engineering Chemistry Research* (2023), DOI `10.1021/acs.iecr.3c02759`.

The main article has now been audited directly. It uses acid-pickled Type-II kerogen prepared from Chang 7 outcrop shale from Tongchuan and therefore remains a **secondary paired dataset**, not a replacement for the locked intact/raw-shale experiment.

### Main-article protocol anchors now recovered

The directly reported experiment defines:

- batch reactor internal volume: `80 cm3`;
- temperature series: `300–700 °C`;
- pressure: `25 MPa`;
- holding time: `2 h`;
- deionized-water : kerogen mass ratio: `1:3`;
- thermocouple accuracy: `±0.5 °C`;
- pressure-sensor accuracy: `±0.05 MPa`;
- initial shale particle size before acid treatment: `120–180 μm`.

The sample preparation is also material to interpretation: carbonate was removed with HCl, silicate/aluminosilicate minerals with HCl + HF, pyrite with HNO3, and free oil by Soxhlet extraction. The resulting solid was treated by the authors as approximately pure kerogen. These differences are why its product yields must not be merged directly with the intact raw-shale values.

### Main-article numerical product anchors now recovered

The directly reported values relevant to the present characterization are:

- peak oil yield: `0.19 g/g TOC` at `380 °C / 25 MPa`;
- light-distillate fraction: `70%` at `500 °C / 25 MPa`;
- asphaltene fraction: `5%` at `500 °C / 25 MPa`;
- methane fraction: `51%` at `600 °C / 25 MPa`;
- hydrogen fraction: `30%` at `700 °C / 25 MPa`;
- gas yield: `0.88 g/g TOC` at `700 °C / 25 MPa`.

The main gas products are reported as CH4, H2, CO2 and C2+ hydrocarbons. Generated oil was characterized by SARA and simulated distillation. The exact extracted rows are stored in `related_pure_kerogen_acs2023.csv` with source location and evidence status.

### What is still not recovered

The ACS article explicitly states that its Supporting Information contains:

- reagents;
- detailed SCW kerogen experimental data;
- detailed gas-production data;
- detailed syngas-component data;
- detailed generated-oil data.

The publisher exposes the SI as a PDF/figshare item, but the complete numerical SI tables have not yet been recovered in this repository audit. In particular, no exact `380 °C` simulated-distillation/boiling-range percentages are entered from the main article. The qualitative trend that gasoline/diesel and other light-distillate fractions rise as temperature is increased toward 500 °C is retained only as `DIRECT_QUALITATIVE` evidence.

Even after the SI is obtained, its values must remain secondary paired pure-kerogen data and cannot replace the primary raw-shale composition.

## 3. Search for same-team / same-sample GC, GC-MS or simulated-distillation data

### Xie et al. 2022 — useful carbon-number evidence, different physical sample

DOI `10.3176/oil.2022.3.02` reports low-maturity organic-rich shale from well F317-181 in the Ordos Basin. The paper reports an Agilent 7890B chromatograph for oil analysis and states that generated oil mainly spans approximately C8-C56. In the 380-450 °C series, fractions below C16 decrease with temperature while fractions above C16 increase.

This is valuable evidence for candidate carbon-number cuts, but it is not the locked Tongchuan outcrop sample: the paper identifies a well sample and reports TOC about 16.25 wt%, versus 15.11 wt% for the primary Tongchuan sample. It must not be used to assign primary-sample lump fractions.

### Xie et al. 2023 Geoenergy Science and Engineering

DOI `10.1016/j.geoen.2023.211553` studies related low-maturity shale in sub/supercritical water, including 360 °C/21 MPa and 400 °C/25 MPa conditions. It is useful for temperature/time trends and future control-case design. Exact identity with the locked Tongchuan physical sample has not been established, so no values are merged into the primary dataset.

### Lu et al. 2026 — same research group / Chang-7 source-rock class, identity not the same by default

DOI `10.1016/j.jaap.2026.107757` includes a Chang-7 Type-II1 source-rock sample and reports a 380 °C oil yield of 234.1 mg/g TOC. This differs materially from the locked primary sample's approximately 352 mg/g TOC, so it cannot be assumed to be the same physical sample or product dataset. It is retained as contextual evidence only until its sample table is matched explicitly.

## 4. Current conclusion

The RSC ESI resolves the primary sample's exact 380 °C SARA and gas-yield observations, including duplicate spread. The ACS paired pure-kerogen main article now adds a verified experimental protocol and a useful `380 °C / 25 MPa` pure-kerogen oil-yield anchor, plus temperature-series quality/gas endpoints.

However, no same-physical-sample 380 °C machine-readable carbon-number distribution, GC/GC-MS oil table, or simulated-distillation table has yet been located for the locked raw-shale product. The remaining critical characterization gap is therefore still the boiling/carbon-number distribution needed to turn the real recovered oil into defensible EOS pseudo-components.
