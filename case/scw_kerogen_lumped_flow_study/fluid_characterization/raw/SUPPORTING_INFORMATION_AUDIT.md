# Supporting-information audit

## Scope

This audit follows the locked research object in `01_RESEARCH_OBJECT.md`: the Tongchuan Chang 7 low-maturity raw shale studied by Zhao et al. at 380 °C, 25 MPa and 4 h. The purpose is to recover primary and clearly separated paired data before any Light/Middle/Heavy lumping.

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

## 2. ACS 2023 paired pure-kerogen study — supplied article bundle audited

Paired study: Zhao et al., *Industrial & Engineering Chemistry Research* (2023), DOI `10.1021/acs.iecr.3c02759`.

The user supplied a Markdown conversion package made from the original 11-page article PDF. Its conversion note states that the article body text, tables, equations, references and embedded raster figures were extracted from the supplied PDF without adding external scientific content. The package includes Figure 1-Figure 11 raster assets, so visible numeric figure labels can be audited directly rather than estimated from axis coordinates.

This article uses acid-pickled Type-II kerogen prepared from Chang 7 outcrop shale from Tongchuan. It therefore remains a **secondary paired dataset**, not a replacement for the locked intact/raw-shale experiment.

### Main-article protocol anchors recovered

The experiment directly defines:

- batch reactor internal volume: `80 cm3`;
- temperature series: `300–700 °C`;
- pressure: `25 MPa`;
- holding time: `2 h`;
- deionized-water : kerogen mass ratio: `1:3`;
- thermocouple accuracy: `±0.5 °C`;
- pressure-sensor accuracy: `±0.05 MPa`;
- initial shale particle size before acid treatment: `120–180 μm`.

Sample preparation is material to interpretation: carbonate was removed with HCl, silicate/aluminosilicate minerals with HCl + HF, pyrite with HNO3, and free oil by Soxhlet extraction. The resulting solid was treated by the authors as approximately pure kerogen. These differences are why its product yields must not be merged directly with intact raw-shale values.

### Table 2 material analysis recovered

The converted article preserves Table 2 exactly. The acid-pickled kerogen row reports:

- C `69.46 wt%`;
- H `6.01 wt%`;
- N `5.58 wt%`;
- S `2.39 wt%`;
- O `11.94 wt%`;
- moisture `2.71 wt%`;
- ash `1.91 wt%`;
- volatile matter `49.16 wt%`;
- fixed carbon `46.22 wt%`.

The original-shale row reports C/H/N/S = `15.84/2.33/2.71/3.80 wt%`. Full rows are stored in `acs2023_table2_material_analysis.csv`.

### 380 °C oil-generation and exact Figure 6 boiling-range anchor

The article reports the pure-kerogen peak oil yield as `0.19 g/g TOC` at `380 °C / 25 MPa`.

Crucially, the supplied Figure 6 raster contains printed values for each simulated-distillation stack. Therefore the 380 °C paired pure-kerogen generated oil is now directly recoverable as:

- gasoline, IBP-180 °C: `0.81 wt%`;
- diesel, 180-350 °C: `23.73 wt%`;
- middle fraction, 350-500 °C: `34.11 wt%`;
- heavy fraction, >500 °C: `41.35 wt%`.

These four values sum to exactly 100%. They are preserved, together with the full 300-500 °C and free-oil Figure 6 series, in `acs2023_figure6_distillation_sara.csv`.

This is a useful secondary-paired constraint for candidate lump boundaries, but it does **not** resolve the primary gap for the locked intact/raw-shale product.

### Figure 6 380 °C SARA source inconsistency

The 380 °C SARA stack in Figure 6 visibly labels:

- saturates `7.21%`;
- aromatics `44.80%` as printed;
- resins `44.45%`;
- asphaltenes `16.78%`.

Those printed labels sum to `113.24%`, so they cannot all represent the same normalized mass-fraction stack. The main text separately states that saturates + aromatics are approximately `38.8%` at the relevant start of the 380-500 °C trend, while resins and asphaltenes are approximately `44.5%` and `16.8%`.

Using the three mutually consistent direct figure labels gives the mass-closure residual:

`100 - 7.21 - 44.45 - 16.78 = 31.56%` aromatics.

Then `7.21 + 31.56 = 38.77%`, matching the rounded main-text `38.8%`. The repository therefore stores `31.56%` only as `DERIVED_MASS_CLOSURE_SOURCE_FIGURE_LABEL_INCONSISTENT`; it is not relabelled as a direct figure measurement. If the ACS SI table is later obtained, it should be used to resolve the apparent figure-label error.

### Additional direct article data recovered

- Figure 6 full simulated-distillation and SARA series for 300, 350, 380, 400, 450 and 500 °C plus free oil;
- Table 4 generated-oil FT-IR indices: at 380 °C, `Xoxid=0.41`, `Xali=0.46`, `Xbrn=0.60`, with the full 380-500 °C series stored in `acs2023_table4_generated_oil_ftir_indices.csv`;
- Table 6 initial/spent-kerogen ultimate analysis and H/C, stored in `acs2023_table6_spent_kerogen_ultimate.csv`;
- Section 3.1 gas-yield anchors: `0.61 g/g TOC` at 650 °C and `0.88 g/g TOC` at 700 °C;
- Section 3.2 carbon-number statements: free-oil C15-C25 `27.73%`; generated-oil C15-C25 spans `17.36-39.51%` over 300-450 °C and falls to `5.57%` at 500 °C;
- Section 3.3 composition anchors: CO2 `16.41%` at 550 °C and `20.62%` at 700 °C; CH4 `41.41%` at 700 °C; the abstract reports the CH4 maximum as about `51%` at 600 °C and H2 about `30%` at 700 °C.

The summary long table is `related_pure_kerogen_acs2023.csv`.

### What remains unrecovered from ACS Supporting Information

The article states that its Supporting Information contains:

- reagents;
- detailed SCW kerogen experimental data;
- detailed gas-production data;
- detailed syngas-component data;
- detailed generated-oil data.

The full SI numerical tables are still not present in the supplied bundle. They remain useful for checking run-by-run values and resolving the Figure 6 SARA inconsistency, but the main article already supplies an exact 380 °C paired pure-kerogen simulated-distillation constraint.

Even after the SI is obtained, its values remain secondary paired pure-kerogen data and cannot replace the primary raw-shale composition.

## 3. Search for same-team / same-sample GC, GC-MS or simulated-distillation data

### Xie et al. 2022 — useful carbon-number evidence, different physical sample

DOI `10.3176/oil.2022.3.02` reports low-maturity organic-rich shale from well F317-181 in the Ordos Basin. The paper reports an Agilent 7890B chromatograph for oil analysis and states that generated oil mainly spans approximately C8-C56. In the 380-450 °C series, fractions below C16 decrease with temperature while fractions above C16 increase.

This is valuable evidence for candidate carbon-number cuts, but it is not the locked Tongchuan outcrop sample: the paper identifies a well sample and reports TOC about 16.25 wt%, versus 15.11 wt% for the primary Tongchuan sample. It must not be used to assign primary-sample lump fractions.

### Xie et al. 2023 Geoenergy Science and Engineering

DOI `10.1016/j.geoen.2023.211553` studies related low-maturity shale in sub/supercritical water, including 360 °C/21 MPa and 400 °C/25 MPa conditions. It is useful for temperature/time trends and future control-case design. Exact identity with the locked Tongchuan physical sample has not been established, so no values are merged into the primary dataset.

### Lu et al. 2026 — same research group / Chang-7 source-rock class, identity not the same by default

DOI `10.1016/j.jaap.2026.107757` includes a Chang-7 Type-II1 source-rock sample and reports a 380 °C oil yield of 234.1 mg/g TOC. This differs materially from the locked primary sample's approximately 352 mg/g TOC, so it cannot be assumed to be the same physical sample or product dataset. It is retained as contextual evidence only until its sample table is matched explicitly.

## 4. Current conclusion

The RSC ESI resolves the primary sample's exact 380 °C SARA and gas-yield observations, including duplicate spread. The supplied ACS article bundle now adds a verified pure-kerogen protocol, exact tabular precursor/oil/spent-kerogen characterization, and—most importantly—an exact visible-label `380 °C / 25 MPa` simulated-distillation profile for the paired pure-kerogen oil.

The primary final-lumping gap is nevertheless unchanged in identity: no **same intact/raw-shale physical sample** 380 °C machine-readable carbon-number distribution or simulated-distillation table has yet been located. The ACS Figure 6 distribution is a strong secondary prior, not a direct replacement.
