# 06 — Scientifically constrained PR parameter baseline

## Purpose

This file defines the first **internally consistent PR parameter baseline** for the zero-flow stage. It is deliberately separated from final calibration: values that come from experimental characterization are fixed as such, values that are literature priors are labelled as priors, and unsupported quantities remain blocked.

The reservoir gate in `05_ZERO_FLOW_PR_CALIBRATION.md` is unchanged. A usable parameter file is not the same thing as a validated model.

## Pure-component / pseudo-component PR inputs

The production PR form uses

`a_i(T) = 0.45724 R^2 Tc_i^2/Pc_i * alpha_i(T)`

and

`b_i = 0.07780 R Tc_i/Pc_i`.

The base component properties are taken from `fluid_characterization/pseudo_component_characterization_380c.csv`; the resulting PR76 kappa, alpha at 653.15 K, `a_c`, `a(T)` and `b` are stored in `pr_parameters/pr_pure_parameters_380c.csv`.

No volume translation is assigned yet. A volume shift is a density parameter and must be fitted independently to density/PVT data; it must never be adjusted through `kij`.

## High-acentric-factor caution

The current production PR path uses the original quadratic PR kappa expression unless the explicit high-omega PR variant is selected. The Diesel, Middle and Heavy pseudo-components have `omega > 0.49`, with Heavy reaching approximately 1.216.

The Peng-Robinson 1978 high-omega extension changes the attraction alpha at 653.15 K by approximately:

- Diesel: 0.12%
- Middle: 1.12%
- Heavy: 3.37%

The comparison is stored in `pr_parameters/pr_alpha_model_sensitivity.csv`.

Therefore PR76 remains the common baseline for immediate regression consistency, but final acceptance of Middle/Heavy requires a PR76-versus-PR78 alpha sensitivity check. A BIP must not be used to compensate for a poor heavy-component alpha function.

## Hydrocarbon–hydrocarbon BIPs

All four oil pseudo-components are contiguous boiling cuts from the same recovered oil. The baseline therefore sets every HC–HC `kij = 0`.

This is not a claim that every real cross interaction is exactly ideal. It is a disciplined starting point: nonzero HC–HC BIPs are introduced only if multicomponent PVT/phase-behavior data demonstrate a systematic error that cannot be explained by pseudo-component characterization.

## H2O–lump initial priors

The standard production PR implementation uses one symmetric `kij`, while the Søreide–Whitson framework distinguishes aqueous and non-aqueous water/hydrocarbon BIPs. Those phase-specific values are therefore used only to define prior information and optimizer seeds; they are not called experimental calibration.

At 653.15 K the Søreide–Whitson pure-water aqueous-form prior evaluated with each pseudo-component's `Tc` and `omega` is approximately:

| lump | SW aqueous prior | SW C5+ non-aqueous prior |
|---|---:|---:|
| Gasoline | 0.0810 | 0.5000 |
| Diesel | 0.0560 | 0.5000 |
| Middle | 0.0966 | 0.5000 |
| Heavy | 0.1509 | 0.5000 |

Because one symmetric PR parameter cannot reproduce two phase-specific Søreide–Whitson values, these numbers define a **model-structure bracket**, not a fitted answer.

For optimizer initialization the baseline chooses the closest available high-T/P proxy evidence:

- `H2O–OIL_GASOLINE`: `kref=0.5000`, `b=0` as a temporary Søreide–Whitson C5+ non-aqueous seed until the n-hexane original rows are recovered.
- `H2O–OIL_DIESEL`: `kref=0.6662345`, `b=-1274.90 K`, reconstructed from the published dodecane PR fits in Teratani et al. This is only an optimizer seed.
- `H2O–OIL_MIDDLE`: `kref=0.2398346`, `b=-477.41 K`, using squalane as the closest saturated equivalent-carbon-number anchor. Aromatic/tetralin priors remain chemistry brackets.
- `H2O–OIL_HEAVY`: the same squalane numbers are retained only as a diagnostic benchmark and `fit_enabled=0`. They are explicitly not the Heavy production BIP.

The prior envelope is recorded in `pr_parameters/water_lump_prior_envelope.csv`, and the full pair matrix in `pr_parameters/pr_binary_matrix_screening.csv`.

## Parameter promotion rules

A water–lump parameter can be promoted from `SCREENING_PRIOR` to `CALIBRATED` only when the corresponding zero-flow dataset passes all of:

1. phase-count topology;
2. both coexistence compositions;
3. density or molar-volume validation;
4. phase-boundary / critical-locus validation;
5. hold-out temperature/pressure validation.

If the optimizer pushes `kij` to a bound, reverses a robust literature temperature trend, or fixes one branch while degrading the other, the result is recorded as a **model-structure failure**, not accepted as a tuned parameter.

## Current baseline matrix at 653.15 K

| pair | kref | b, K | status |
|---|---:|---:|---|
| H2O–Gasoline | 0.5000 | 0 | screening prior |
| H2O–Diesel | 0.66623 | -1274.90 | dodecane proxy prior |
| H2O–Middle | 0.23983 | -477.41 | squalane EACN prior |
| H2O–Heavy | 0.23983 | -477.41 | benchmark only; fitting forbidden |
| any oil-lump pair | 0 | 0 | locked screening baseline |

These values are scientifically traceable **starting parameters**, not reservoir-ready calibration results.
