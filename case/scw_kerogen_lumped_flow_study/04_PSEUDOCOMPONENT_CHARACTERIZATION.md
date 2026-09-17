# 04 — Pseudo-component characterization at 380 °C

## Decision

The first EOS-ready characterization is built from the **measured 380 °C recovered-oil boiling-range fractions plus the measured/digitized 380 °C carbon-number profile**. It is not built by choosing `nC4`, `nC10`, or squalane and copying their pure-component properties.

The resulting table is suitable as an **initial PVT/EOS screening characterization**. It is not a claim that every value was measured directly. Each row keeps the experimental inputs separate from the petroleum-characterization estimates.

Machine-readable results are in:

`fluid_characterization/pseudo_component_characterization_380c.csv`

The SCN reconstruction used to form the average properties is in:

`fluid_characterization/raw/figure5_380c_scn_characterization_basis.csv`

## Characterization table

All compositions below are on the **four-lump recovered generated-oil basis**. They do not yet include H2O or the separately recovered gas products.

| pseudo-component | boiling range | mass fraction | mole fraction | MW, g/mol | average Tb, °C | SG 60/60 | Tc, K | Pc, MPa | omega | Vc, cm3/mol |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `OIL_GASOLINE` | IBP–180 °C | 0.008100 | 0.036330 | 82.392 | 68.991 | 0.691443 | 515.231 | 3.192179 | 0.269328 | 358.490 |
| `OIL_DIESEL` | 180–350 °C | 0.237300 | 0.406770 | 215.583 | 278.493 | 0.836212 | 736.423 | 1.751155 | 0.583473 | 824.073 |
| `OIL_MIDDLE` | 350–500 °C | 0.341100 | 0.325421 | 387.349 | 430.304 | 0.896458 | 870.437 | 1.138075 | 0.896554 | 1318.853 |
| `OIL_HEAVY` | >500 °C | 0.413500 | 0.231479 | 660.132 | 563.696 | 0.942605 | 982.864 | 0.808714 | 1.215676 | 1793.890 |

Mass fractions sum to `1.000000` and the MW-derived mole fractions sum to `1.000000`.

## Experimental inputs

### 1. Lump mass fractions

The four mass fractions are direct Figure 6 labels from Zhao et al., *Industrial & Engineering Chemistry Research* 2023, DOI `10.1021/acs.iecr.3c02759`, for generated oil at 380 °C and 25 MPa:

- IBP–180 °C: `0.81 wt%`;
- 180–350 °C: `23.73 wt%`;
- 350–500 °C: `34.11 wt%`;
- >500 °C: `41.35 wt%`.

These are `SECONDARY_PAIRED` because the experiment uses acid-pickled Type-II kerogen rather than the locked intact/raw-shale product.

### 2. Carbon-number structure inside each lump

Figure 5 of the same paper gives the 380 °C carbon-number distribution. The raster was digitized only to obtain **relative SCN weights inside each already-fixed boiling-range lump**. Figure 5 is not used to overwrite the Figure 6 mass fractions.

The SCN-to-lump mapping follows the Katz–Firoozabadi generalized average boiling-point table:

- `OIL_GASOLINE`: `<C7` proxy through C10; C10 average Tb is below 180 °C while C11 is above it;
- `OIL_DIESEL`: C11–C20; C20 average Tb is about 338 °C while C21 is about 351 °C;
- `OIL_MIDDLE`: C21–C37; C37 average Tb is about 500 °C;
- `OIL_HEAVY`: C38+; the visible Figure 5 tail extends to approximately C74.

Because the 380 °C Figure 5 panel contains only very small visible sub-C11 bars, `OIL_GASOLINE` has the largest relative digitization uncertainty even though it contributes only 0.81 wt% of the recovered oil.

## Petroleum-characterization methods

### Molecular weight

For SCN `I`, the Pedersen-style generalized petroleum relation is used:

`MW_I = 14.0269 I - 4`

The lump average MW is obtained from the digitized **mass** shares by the harmonic mass-to-mole conversion:

`MW_lump = 1 / sum(w_I / MW_I)`.

The four-lump mole fraction then follows from:

`z_j = (W_j / MW_j) / sum_k(W_k / MW_k)`.

This means the reported mole fractions are derived characterization values, not directly measured compositions.

### Normal boiling point and specific gravity

For C6–C45, SCN average boiling points and specific gravities use the generalized Katz–Firoozabadi petroleum-fraction table. That table was developed from multiple condensate/crude-oil systems and is used here as a petroleum characterization basis, not as proof that the products are normal alkanes.

For the visible C46+ tail:

- SCN boiling point is extended with the Twu normal-paraffin reference relation;
- SCN specific gravity is log-extrapolated from the C30–C45 generalized petroleum trend;
- lump Tb is the mass-weighted SCN mean;
- lump SG uses ideal-volume mass mixing: `SG_lump = 1 / sum(w_I / SG_I)`.

The C46+ extrapolation is one reason `OIL_HEAVY` is explicitly provisional.

### Critical properties and acentric factor

`Tc`, `Pc`, and `Vc` are calculated using the internally consistent Twu (1984) petroleum-fraction correlation from lump normal boiling point and specific gravity:

C. H. Twu, “An internally consistent correlation for predicting the critical properties and molecular weights of petroleum and coal-tar liquids,” *Fluid Phase Equilibria* 16 (1984) 137–150, DOI `10.1016/0378-3812(84)85027-X`.

Twu's reported applicability covers boiling points and specific gravities beyond the four current lump averages, so these estimates are within the published input range. The acentric factor is then calculated with the Edmister correlation from `Tb`, `Tc`, and `Pc`.

These correlations produce internally usable cubic-EOS starting properties. They do **not** replace calibration against high-temperature/high-pressure VLE/LLE/PVT data.

## Heavy treatment

`OIL_HEAVY` must be treated differently from the lighter fractions.

The central characterization in the table uses the visible C38–C74 Figure 5 tail and returns approximately:

- `MW = 660.1 g/mol`;
- `Tb = 563.7 °C`;
- `SG = 0.9426`;
- `Tc = 982.9 K`;
- `Pc = 0.809 MPa`;
- `omega = 1.216`;
- `Vc = 1794 cm3/mol`.

This row is marked `PROVISIONAL_HEAVY_TAIL_LOWER_BOUND_LIKE`, not as an exact measured heavy-pseudocomponent property set. The reason is physical: Figure 6 assigns 41.35 wt% of the recovered oil to >500 °C material, while SARA shows substantial resin/asphaltene character. The GC-visible carbon-number tail can under-represent the least volatile and most polar material. Missing material would generally push the effective Heavy MW/Tb/SG characterization upward or otherwise change its PNA/associating character.

Therefore:

1. do not tune away this uncertainty by forcing `OIL_HEAVY` to a convenient pure hydrocarbon;
2. keep Heavy MW/SG/Tb and CPA association behavior as high-priority validation targets;
3. update the Heavy row when measured heavy-cut MW, density, TBP/simdist residue information, or high-T/P phase-equilibrium data become available;
4. preserve the experimental >500 °C cut even if the property values are later revised.

## Squalane policy: benchmark only

Squalane remains useful as a **heavy saturated-hydrocarbon benchmark**, not as the identity of `OIL_HEAVY`.

Reference values for the benchmark are stored in `fluid_characterization/squalane_benchmark.csv`. The benchmark has approximately:

- formula `C30H62` and MW `422.8133 g/mol` (NIST Chemistry WebBook);
- normal boiling point `350 °C` and density `0.8115` at 15 °C (PubChem/PAC compilation);
- measured critical temperature `822 K` and critical pressure `0.700 MPa` by Nikitin and Popov, *Fluid Phase Equilibria* 237 (2005) 16–20, DOI `10.1016/j.fluid.2005.06.026` (NIST ThermoML).

Those values are very different from the current >500 °C experimental Heavy characterization, especially in MW, boiling point, and density. Squalane can therefore be retained for code regressions, sensitivity studies, and a nonpolar-heavy reference, but it must not silently populate the production `OIL_HEAVY` row.

## Data-status policy

The table intentionally distinguishes three evidence layers:

- **direct experiment:** Figure 6 lump mass fractions and Figure 5 carbon-number shape;
- **petroleum characterization:** SCN MW/Tb/SG mapping and Twu/Edmister pseudo-critical estimates;
- **unresolved calibration:** primary raw-shale fraction weights, Heavy residue characterization, H2O–pseudo-component interaction parameters, CPA association parameters, and high-T/P viscosity/phase-equilibrium validation.

For the next thermodynamic step, the four rows may be used as the **initial PR/SW/CPA screening property set**, but they must retain their provenance/status fields. Final model claims require phase-behavior calibration rather than treating these correlations as experimental critical-property measurements.

## References used for characterization

1. Zhao et al., *Industrial & Engineering Chemistry Research* 62 (2023) 17343–17353, DOI `10.1021/acs.iecr.3c02759` — Figure 5 carbon-number distribution and Figure 6 simulated distillation.
2. Katz and Firoozabadi, “Predicting Phase Behavior of Condensate/Crude-Oil Systems Using Methane Interaction Coefficients,” *Journal of Petroleum Technology* 30 (1978), DOI `10.2118/6721-PA` — generalized SCN petroleum properties.
3. Pedersen petroleum-plus-fraction characterization — generalized SCN MW relation used for the Figure 5 reconstruction.
4. Twu, *Fluid Phase Equilibria* 16 (1984) 137–150, DOI `10.1016/0378-3812(84)85027-X` — internally consistent petroleum-fraction critical-property correlation.
5. Edmister acentric-factor correlation — applied after Twu `Tc/Pc` characterization.
6. Nikitin and Popov, *Fluid Phase Equilibria* 237 (2005) 16–20, DOI `10.1016/j.fluid.2005.06.026` — measured squalane critical properties; benchmark only.
