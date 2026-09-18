# Supporting-information audit

## Scope

This audit follows the locked research object used by this case: the Tongchuan Chang 7 low-maturity raw shale studied by Zhao et al. at 380 °C, 25 MPa and 4 h. The purpose is to recover primary supplementary data before any Light/Middle/Heavy lumping.

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

## 2. ACS 2023 Supporting Information — recovered and audited

Paired study: Zhao et al., *Industrial & Engineering Chemistry Research* (2023), DOI `10.1021/acs.iecr.3c02759`.

The study uses acid-pickled Type-II kerogen rather than the locked intact raw-shale sample. Its values are therefore `SECONDARY_PAIRED`: they may constrain candidate boiling-range cuts and sensitivity ranges, but may not replace the Tongchuan raw-shale mass fractions required by M1.

### 2.1 File provenance

The public ACS Figshare records are:

- collection: `6876126`;
- journal-contribution record: `24291372`;
- file ID: `42637827`;
- recovered filename: `ie3c02759_si_001.pdf`;
- SHA-256: `ed2f1396f44fac84389bf4f165e1b89603dfb6d8585a5b3697c397767975f117`;
- PDF pages: 6.

The source PDF is external literature and is not committed to this repository. Only traceable numeric extractions and the source hash are retained.

### 2.2 380 °C experimental state from SI Table S2

The nominal state is 380 °C and 25 MPa with a common reaction time of 2 h and designed water/kerogen mass ratio of 3.

The two 380 °C experiments report:

| replicate | measured T (°C) | measured P (MPa) | kerogen (g) | water (g) | measured water/kerogen ratio |
|---|---:|---:|---:|---:|---:|
| 1 | 379.8 | 24.70 | 6.40 | 19.29 | 3.01 |
| 2 | 382.6 | 24.92 | 6.40 | 19.30 | 3.02 |

These conditions differ materially from the locked RSC raw-shale experiment, which used intact shale, approximately 1:1 water/shale mass ratio, and 4 h reaction time.

### 2.3 Exact 380 °C generated-oil data from SI Table S5

The source reports two generated-oil replicates and an average:

| quantity | replicate 1 | replicate 2 | source average |
|---|---:|---:|---:|
| recovered oil mass (g) | 0.84 | 0.82 | 0.83 |
| oil yield (g/g kerogen) | 0.13 | 0.13 | 0.13 |
| oil yield (g/g TOC) | 0.19 | 0.18 | 0.19 |
| saturates (wt%) | 6.42 | 7.98 | 7.20 |
| aromatics (wt%) | 30.82 | 32.30 | 31.56 |
| resins (wt%) | 45.35 | 43.56 | 44.45 |
| asphaltene-equivalent column (wt%) | 17.41 | 16.16 | 16.78 |

Table S5 literally labels its fourth SARA-like column `Aliphatene`. The mapping above records it as an asphaltene-equivalent column because the four columns close to 100%, and the main article's rounded statement of about 5% asphaltene at 500 °C corresponds to the Table S5 average of 5.51% in the same column. The original header spelling is preserved in `related_pure_kerogen_acs2023.csv` notes.

### 2.4 Exact 380 °C simulated-distillation fractions from SI Table S5

| fraction | replicate 1 (wt%) | replicate 2 (wt%) | source average (wt%) |
|---|---:|---:|---:|
| Gasoline | 0.79 | 0.83 | 0.81 |
| Diesel | 22.68 | 24.78 | 23.73 |
| Distillate | 33.97 | 34.25 | 34.11 |
| Heavy Oil | 42.56 | 40.14 | 41.35 |
| **sum** | **100.00** | **100.00** | **100.00** |

The source-average non-heavy sum is therefore `0.81 + 23.73 + 34.11 = 58.65 wt%`. This is a derived check, not a separate source field.

The SI table gives the four product-fraction labels and percentages but does not itself state the boiling-point or carbon-number boundaries associated with those labels. Therefore these fractions are not yet a defensible direct mapping to the case's candidate `C6-C14 / C15-C20 / C21+` lumps.

### 2.5 SimDist cut-definition audit and candidate carbon-number mapping

The recovered ACS SI Table S5 reports the labels `Gasoline / Diesel / Distillate / Heavy Oil` and their exact fractions, but the SI itself does **not** print numerical boiling-point cut temperatures. The accessible ACS article page also does not expose the paid full-text method/figure containing a direct numeric definition. Therefore the cut temperatures below are not recorded as `DIRECT_ACS_NUMERIC`.

A high-confidence convention can nevertheless be established from an evidence chain:

1. a later open-access paper from the same research group (Tian Xie, Hui Jin and Qiuyang Zhao among the authors; DOI `10.3176/oil.2025.2.02`) again reports the same four oil fractions and explicitly describes its figure in terms of `IBP` (initial boiling point) and `FBP` (final boiling point);
2. shale-oil simulated-distillation literature explicitly defines the corresponding bins as gasoline below 180 °C, diesel 180-350 °C, the 350-500 °C distillate/VGO interval, and heavy oil above 500 °C;
3. an independent heavy-oil distillation table using the exact `Gasoline / Diesel / Distillate / Heavy oil` labels gives `IBP-180 / 180-350 / 350-500 / >500 °C`.

For repository use this is therefore classified as `HIGH_CONFIDENCE_INFERRED_CONVENTION`, not a direct ACS Table S5 field:

| ACS label | candidate boiling range | n-alkane-equivalent carbon range | basis |
|---|---:|---|---|
| Gasoline | IBP-180 °C | IBP to approximately C10 | n-C10 boils at about 174 °C and n-C11 at about 195 °C |
| Diesel | 180-350 °C | approximately C11-C20 | n-C20 boils at about 343 °C and n-C21 at about 357 °C |
| Distillate | 350-500 °C | approximately C21-C36 | 500 °C lies approximately between n-C36 and n-C37 |
| Heavy Oil | >500 °C | approximately C37+ | n-C37 is just above 500 °C |

These carbon numbers are **normal-alkane-equivalent boiling indices**, not chemical identities. This distinction matters for the recovered oil because its SARA data show large aromatic/resin/asphaltene fractions.

The previous case proposal `C6-C14 / C15-C20 / C21+` is not fully observable from the recovered ACS SimDist bins. The `180-350 °C` Diesel interval spans approximately C11-C20 and crosses the proposed C14/C15 boundary, so it cannot be split into experimental C6-C14 and C15-C20 mass fractions without carbon-number-resolved data. By contrast, the previous broad `C21+` boundary is approximately compatible with the 350 °C cut and can be formed as `Distillate + Heavy Oil` for **secondary sensitivity analysis only**.

A SimDist-aligned three-lump candidate is therefore:

| candidate lump | boiling range | approximate carbon range | replicate 1 (wt%) | replicate 2 (wt%) | mean (wt%) | observed range | half-range |
|---|---:|---|---:|---:|---:|---:|---:|
| Light | IBP-180 °C | <=C10 equivalent | 0.79 | 0.83 | 0.81 | 0.79-0.83 | +/-0.02 |
| Middle | 180-350 °C | C11-C20 equivalent | 22.68 | 24.78 | 23.73 | 22.68-24.78 | +/-1.05 |
| Heavy | >350 °C | C21+ equivalent | 76.53 | 74.39 | 75.46 | 74.39-76.53 | +/-1.07 |

For the original four source bins, the duplicate half-ranges are Gasoline +/-0.02, Diesel +/-1.05, Distillate +/-0.14 and Heavy Oil +/-1.21 wt%. With only two replicates, the repository uses raw ranges/half-ranges rather than treating a sample standard deviation as a well-characterized experimental uncertainty.

The machine-readable mapping and duplicate statistics are stored in `acs2023_simdist_lumping_candidates.csv`.

### 2.5 What the recovered ACS SI changes

The ACS access gap is closed: exact 380 °C generated-oil and distillation data are now available and retained in `related_pure_kerogen_acs2023.csv`.

It does **not** close the primary M1 mass-fraction gate because:

1. the sample is acid-pickled kerogen rather than the locked Tongchuan intact raw shale;
2. the reaction time and water/feed ratio differ;
3. Table S5's distillation category boundaries are not stated in the SI itself;
4. M1 requires experimental mass fractions for the final lumps of the locked product, not merely a related kerogen product.

## 3. Search for same-team / same-sample GC, GC-MS or simulated-distillation data

### Xie et al. 2022 — useful carbon-number evidence, different physical sample

DOI `10.3176/oil.2022.3.02` reports low-maturity organic-rich shale from well F317-181 in the Ordos Basin. The paper reports an Agilent 7890B chromatograph for oil analysis and states that generated oil mainly spans approximately C8-C56. In the 380-450 °C series, fractions below C16 decrease with temperature while fractions above C16 increase.

This is valuable evidence for candidate carbon-number cuts, but it is not the locked Tongchuan outcrop sample: the paper identifies a well sample and reports TOC about 16.25 wt%, versus 15.11 wt% for the primary Tongchuan sample. It must not be used to assign primary-sample lump fractions.

### Xie et al. 2023 Geoenergy Science and Engineering

DOI `10.1016/j.geoen.2023.211553` studies related low-maturity shale in sub/supercritical water, including 360 °C/21 MPa and 400 °C/25 MPa conditions. It is useful for temperature/time trends and future control-case design. Exact identity with the locked Tongchuan physical sample has not been established, so no values are merged into the primary dataset.

### Lu et al. 2026 — same research group / Chang-7 source-rock class, identity not the same by default

DOI `10.1016/j.jaap.2026.107757` includes a Chang-7 Type-II1 source-rock sample and reports a 380 °C oil yield of 234.1 mg/g TOC. This differs materially from the locked primary sample's approximately 352 mg/g TOC, so it cannot be assumed to be the same physical sample or product dataset. It is retained as contextual evidence only until its sample table is matched explicitly.

### 2026-09-18 targeted literature re-check

Targeted searches using the locked Tongchuan identifiers (`Tongchuan`, `TOC 15.11 wt%`, Zhao/Dong/Xie/Jin/Guo authorship, the exact RSC DOI and the reported 352.1 mg/g TOC oil-yield anchor) did not locate an additional publication that exposes a machine-readable 380 °C carbon-number distribution or simulated-distillation table for the same physical raw-shale sample.

This is evidence of a search attempt, not proof that unpublished or inaccessible data do not exist. The acceptance rule remains: a different Chang-7 sample cannot be promoted to primary evidence by similarity alone.

## 4. Current conclusion

The RSC ESI resolves the primary sample's exact 380 °C SARA and gas-yield observations, including duplicate spread. The recovered ACS SI now adds exact paired-kerogen 380 °C SARA and simulated-distillation constraints.

However, no same-physical-sample 380 °C machine-readable carbon-number distribution, GC/GC-MS oil table, or simulated-distillation table has yet been located for the locked Tongchuan raw-shale product.

Therefore the remaining critical characterization gap is still the boiling/carbon-number distribution needed to turn the real recovered oil into defensible EOS pseudo-components. The ACS dataset narrows plausible sensitivity space but does not make M1 pass.
