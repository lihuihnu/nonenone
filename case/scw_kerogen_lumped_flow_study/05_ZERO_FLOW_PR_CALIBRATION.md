# 05 — Zero-flow PR calibration before reservoir flow

## Hard decision gate

**Reservoir flow is frozen until all four H2O–pseudo-component binary PR calibrations pass their zero-flow validation gates.**

The calibrated object is not “a `kij` that makes flash converge.” The calibrated object is a binary phase-behavior model that must reproduce, on independent experimental information:

1. phase count / topology;
2. both coexisting-phase compositions;
3. phase or homogeneous-mixture density/volume information relevant to the two branches;
4. phase-boundary / critical-locus position.

A candidate `kij(T)` that improves nonlinear convergence but fails any required physical observable is rejected.

## Model form

The first PR calibration uses the existing production Peng–Robinson EOS and production flash implementation. The default temperature form is

`kij(T) = k_ref + b * (1/T - 1/T_ref)`, with `T_ref = 653.15 K`.

This form is fitted only when at least two distinct temperatures have admissible raw experimental information. A constant `kij` is allowed for a single-temperature diagnostic but cannot be promoted as the production temperature law.

`lij` or a non-classical mixing rule is **not** silently introduced to improve a fit. If classical PR plus one `kij(T)` fails the validation gate, the failure is recorded as model-structure evidence.

## Experimental-proxy policy

The pseudo-components are petroleum fractions, not pure compounds, so no literal experimental binary named `H2O–OIL_MIDDLE` exists. Calibration therefore uses experimentally measured proxy families selected by the lump's boiling range, equivalent carbon number and chemical character.

A proxy is admissible only as far as its experimental coverage supports the target lump. A fitted proxy `kij` is never copied blindly into the pseudo-component: the production pseudo-component fit must be reconstructed with the repository's own `Tc/Pc/omega` and evaluated against the proxy phase data.

The source inventory and admissibility status are machine-readable in `binary_pr_calibration/source_manifest.csv`.

## Binary-by-binary plan and current status

### H2O–OIL_GASOLINE

Primary high-T/P anchor: water + n-hexane.

Tian, Michelberger and Franck (1991), DOI `10.1016/S0021-9614(05)80063-1`, experimentally determined the two-phase surface for water + n-hexane over approximately `550–700 K` and `20–240 MPa`, including the water-critical region. The paper also reports homogeneous-mixture molar volumes and a critical-curve minimum near `628 K / 31 MPa`.

This is excellent coverage for phase count, phase-boundary location and volumetric behavior around the `653 K / 25 MPa` target. However, the exact machine-readable table rows needed for branch-by-branch regression have not yet been committed in this case. Therefore:

- phase-boundary source: **located**;
- volume/density source: **located**;
- exact branch-composition table for production regression: **not yet ingested**;
- production `kij(T)`: **not accepted yet**.

### H2O–OIL_DIESEL

Primary anchor: water + dodecane; decane data are a secondary lighter anchor.

Stevenson et al. (1994), DOI `10.1016/0378-3812(94)87016-0`, measured VLE/LLE compositions and critical points for water + dodecane at `600–660 K` and pressures up to about `31 MPa`, directly overlapping the target temperature/pressure neighborhood.

Teratani et al. (2017), DOI `10.1627/jpi.60.26`, published PR fits for water + dodecane at `603.6 K` and `633.0 K`, and for water + decane at several VLE/LLE temperatures. Those published BIPs are retained only as **literature priors / regression checks**, because their pure-component characterization and mixing implementation are not identical to this repository.

Current status:

- coexistence composition + critical-locus source: **located**;
- exact raw dodecane tables: **not yet ingested into this case**;
- direct binary phase-density source near target: **not yet adequate**;
- production `kij(T)`: **not accepted yet**.

### H2O–OIL_MIDDLE

Hydrocarbon-number anchor: water + squalane, supplemented by aromatic/naphthenic proxy evidence such as 1-methylnaphthalene and tetralin.

Squalane is near the middle lump's equivalent carbon number, and Stevenson et al. provide high-T/P VLE/LLE information. The repository already contains exact water–squalane coexistence rows at `637.2 K` and `653.2 K` plus hold-out points.

However, the existing zero-flow PR regression is an explicit **failure of the full gate**: it can improve the hydrocarbon-rich composition but does not reproduce both branches within experimental tolerances, and the associated density treatment remains outside the reference uncertainty. Therefore that fit is a model-structure diagnostic, not an accepted `OIL_MIDDLE` BIP.

Aromatic/naphthenic PR fits from Teratani et al. are useful as chemistry brackets, because the generated oil is not a pure saturated paraffin fraction.

Current status: **partial proxy data; PR structure/chemistry validation still open**.

### H2O–OIL_HEAVY

This is the strictest gate.

`OIL_HEAVY` is an experimental `>500 °C` fraction with substantial resin/asphaltene character. Squalane is retained only as a saturated-heavy benchmark and is **not an admissible identity substitute** for this pseudo-component.

Useful evidence includes:

- Stevenson et al. water–squalane data as a nonpolar lower-complexity benchmark;
- Sato et al. (2018), DOI `10.1627/jpi.61.256`, measured water + atmospheric-residue VLE at `603–643 K` and `2.0–10.2 MPa`, providing chemically closer heavy-oil evidence, but not the 25 MPa target LLE/density coverage.

The Heavy `kij(T)` is blocked until a defensible data envelope includes heavy-rich and water-rich equilibrium information plus phase-boundary and density/volume validation. **No Heavy BIP may be changed merely to obtain a two-phase flash or improve convergence.**

Current status: **BLOCKED — no adequate target-window experimental proxy envelope yet**.

## Regression and validation protocol

For each binary proxy dataset:

1. preserve original rows and reported experimental uncertainties; do not invent standard deviations;
2. split calibration/validation by temperature and/or pressure bands, not by random individual rows;
3. fit `k_ref` and `b` only on calibration rows using production PR fugacity/flash paths;
4. evaluate every held-out state for phase count before composition error is scored;
5. compare both equilibrium branches, not only the hydrocarbon-rich branch;
6. compare available experimental density/molar-volume information without changing `kij` to compensate for a pure-component density error;
7. compare the predicted phase boundary / critical locus against experimental boundary data;
8. accept only if all mandatory gates pass with pre-declared tolerances.

A volume-translation parameter, if used, is fitted independently to density data and must not change fugacity equilibrium. `kij` is not allowed to absorb a density-model error.

## Required acceptance matrix

Each H2O–lump binary must reach `PASS` for:

| Gate | Requirement |
|---|---|
| Phase count | Correct 1/2/3-phase classification over calibration and hold-out boundary states |
| Two-phase compositions | Both coexisting branches reproduced within reported uncertainty or an explicitly justified engineering tolerance |
| Phase density / volume | Independent experimental density or molar-volume target reproduced; pure-component-only density is insufficient for final acceptance |
| Phase boundary | Boundary/critical-locus position reproduced without branch switching or spurious islands |
| Hold-out behavior | No acceptance based solely on training rows |
| Parameter physicality | Smooth `kij(T)`; no discontinuous point-by-point tuning |

The current machine-readable gate is `binary_pr_calibration/reservoir_entry_gate.csv`.

## Current reservoir-entry decision

**BLOCKED.**

The repository must not create or promote the final reservoir flow case from this pseudo-component set until all four binary rows pass. The existing H2O–squalane fit is specifically *not* a Heavy pass and is not sufficient for Middle either.

## Literature anchors

- Tian, Y.; Michelberger, T.; Franck, E. U. (1991), *J. Chem. Thermodynamics* 23, 105–112, DOI `10.1016/S0021-9614(05)80063-1` — water + n-hexane phase boundary and molar-volume data.
- Stevenson, R. L.; LaBracio, D. S.; Beaton, T. A.; Thies, M. C. (1994), *Fluid Phase Equilibria* 93, 317–336, DOI `10.1016/0378-3812(94)87016-0` — water + dodecane and water + squalane VLE/LLE compositions and critical phenomena.
- Teratani, S.; Ota, M.; Sato, Y.; Inomata, H. (2017), *J. Jpn. Petrol. Inst.* 60, 26–33, DOI `10.1627/jpi.60.26` — PR correlation study and literature BIP priors for several water–hydrocarbon systems.
- Sato, S. et al. (2018), *J. Jpn. Petrol. Inst.* 61, 256–262, DOI `10.1627/jpi.61.256` — water + atmospheric-residue VLE and heavy-oil structural characterization.
