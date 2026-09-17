# 08 — Independent density and viscosity validation sets

## Decision

Density and viscosity validation are separated from phase-equilibrium calibration and from each other.

Neither validation set may be used to fit `kij`. The validation CSV rows are also not allowed to fit volume translation or viscosity-correlation parameters. If a density shift or transport parameter requires fitting, that fit must use a different calibration subset/source and then be evaluated against the files below without retuning.

## Density validation set

Machine-readable rows:

`property_validation/density/density_validation.csv`

The current accepted rows contain two evidence layers:

1. **H2O formulation verification** — official IAPWS-IF97 Region-3 verification states. These test implementation correctness of the aqueous density closure.
2. **Squalane independent hold-out benchmark** — only the rows that the historical squalane calibration already labelled `validation` are reused. Historical training rows are excluded so they cannot masquerade as independent validation.

Squalane remains a saturated-heavy benchmark for Middle/Heavy, not identity data for either pseudo-component.

Additional independent density sources have been located but not yet promoted into pointwise rows:

- Span–Wagner n-hexane EOS for the Gasoline/light bracket;
- Wu et al. n-C16/n-C18/n-C20 high-pressure density measurements to 523 K and 265 MPa;
- Romeo–Lemmon n-C16/n-C22 reference EOS.

The source inventory and row-ingestion state are in `property_validation/density/source_manifest.csv`.

### Density acceptance philosophy

For water, official verification values are a numerical implementation gate.

For pseudo-components, matching a proxy is necessary but not sufficient. A proxy row can reject an obviously poor density model, but it cannot by itself establish that a pseudo-component is experimentally validated.

Volume translation is explicitly outside this validation fit. A proposed shift must be frozen before these rows are evaluated.

## Viscosity validation set

Machine-readable rows:

`property_validation/viscosity/viscosity_validation.csv`

The transport validation is deliberately performed in two stages.

### Stage 1 — intrinsic transport closure

Use the **reference density in the validation row** and call the density-based viscosity closure. This isolates the transport model:

`mu_model = f(T, rho_reference, composition)`

rather than contaminating the result with EOS density error.

The current rows include:

- official IAPWS-2008 water verification points;
- n-hexane reference-correlation computer-verification points for the Gasoline bracket;
- n-undecane reference-correlation verification points as a lower-EACN Diesel bracket;
- the independent squalane hold-out rows as a saturated-heavy benchmark.

The n-hexane reference correlation ends at 600 K, so it does **not** validate the 653 K target. The dense n-undecane correlation is validated over a more limited liquid range; its 635 K dense point is retained only as an extrapolation diagnostic.

### Stage 2 — coupled EOS + transport

Only after Stage 1 is acceptable should the model replace `rho_reference` with EOS-predicted density. This second test measures the accumulated density + transport error and must be reported separately.

A model that fails Stage 1 cannot be rescued by tuning the EOS density.

## Long-chain transport data still to ingest

Baled et al. provide n-C16/n-C18/n-C20 viscosities over 304–534 K and 3–243 MPa. NIST ThermoML exposes the raw datasets and they are the preferred next addition for Diesel/Middle.

Pimentel-Rodas et al. provide simultaneous density and viscosity for pentane/octane/nonane/decane/dodecane to 30 MPa and are useful as a lower-temperature consistency set.

These sources are listed as `ROWS_PENDING`; their ranges are not converted into synthetic points.

## Heavy policy

Heavy has no independent pass yet.

The squalane rows are deliberately marked benchmark-only. They cannot certify the real >500 °C fraction, whose MW, aromatic/resin/asphaltene content and polarity differ substantially.

Heavy density and viscosity remain blocked until a heavy-cut/residue dataset or a defensible reconstructed-oil experiment is available. Neither PR nor CPA may use a density shift, critical volume, LBC parameter, CPA association term or BIP to hide this missing validation.

## Gates

- density gate: `property_validation/density/density_gate.csv`
- viscosity gate: `property_validation/viscosity/viscosity_gate.csv`

Current overall status:

- `DENSITY_GATE_BLOCKED`
- `VISCOSITY_GATE_BLOCKED`

This is intentional and keeps transport/property uncertainty separate from phase-equilibrium calibration.
