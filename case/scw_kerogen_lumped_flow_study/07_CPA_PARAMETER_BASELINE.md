# 07 — Independent CPA parameter baseline

## Decision

CPA is treated as a **separate thermodynamic model**, not as PR with an association correction pasted on top. The model owns an independent parameter set:

- SRK-family physical parameters `a0`, `b`, and `c1`;
- association scheme and site counts;
- self-association energy `epsilon` and volume `beta`;
- optional cross-association parameters;
- CPA-specific binary interaction parameters `kij`.

PR `kij(T)`, PR `a/b`, and PR alpha-model choices are not imported.

## Water parameter set

The baseline water model is the 4C CPA parameter set already exercised by the repository's public CO2–DME–H2O VLLE benchmark:

- `a0 = 0.12277 Pa m6/mol2`;
- `b = 1.4515e-5 m3/mol`;
- `c1 = 0.67359`;
- two donor + two acceptor sites;
- `epsilon = 16655 J/mol`;
- `beta = 0.0692`;
- SRK physical term;
- simplified CPA radial distribution in the current implementation.

This is the one part of the present CPA parameter system that is already tied to a repository-level CPA validation case.

## Petroleum pseudo-component physical parameters

For baseline v1, all four oil pseudo-components use one internally consistent mapping from the experiment-driven `Tc/Pc/omega` characterization into the SRK physical term:

`a0 = 0.42748 R^2 Tc^2/Pc`

`b = 0.08664 R Tc/Pc`

`c1 = 0.480 + 1.574 omega - 0.176 omega^2`.

This gives:

| component | a0, Pa m6/mol2 | b, m3/mol | c1 |
|---|---:|---:|---:|
| Gasoline | 2.45754 | 1.16270e-4 | 0.89116 |
| Diesel | 9.15196 | 3.02939e-4 | 1.33847 |
| Middle | 19.67379 | 5.50959e-4 | 1.74971 |
| Heavy | 35.30011 | 8.75490e-4 | 2.13337 |

These are **screening physical parameters**, not vapor-pressure/density-fitted CPA pure parameters. Yan, Kontogeorgis and Stenby (2009) explicitly developed CPA narrow-cut correlations for ill-defined C7+ reservoir fractions; that methodology is the preferred future upgrade for Diesel/Middle/Heavy once its equations and validation basis are fully reproduced in this repository.

## Association policy for oil lumps

The baseline does not assign self-association to any oil lump.

This is deliberate. CPA reduces to the SRK physical term for non-self-associating species, and water–alkane systems are routinely modeled this way. Aromatics may require water cross-association/solvation, but an experimental mixed petroleum cut cannot be given a benzene-like site count or cross-association beta without a composition model.

Therefore:

- Gasoline: no association;
- Diesel: no association;
- Middle: no self-association; aromatic solvation remains an explicit open model option;
- Heavy: no self-association in baseline v1; polar association is **unresolved**, not zero by physical assertion.

The Heavy SARA evidence is not sufficient to invent donor/acceptor counts, `epsilon`, or `beta`. If CPA eventually requires an associating heavy representation, the preferred route is to obtain PNA/SARA/heteroatom and phase data and consider splitting Heavy into nonpolar and polar pseudo-components, rather than forcing one mixed >500 C lump into an arbitrary association scheme.

## CPA-specific water–hydrocarbon BIP priors

Oliveira, Coutinho and Queimada (2007) used standard CPA with one temperature-independent physical `kij` for n-alkane + water systems. Their directly tabulated values include:

- n-hexane: `kij = 0.044`;
- n-decane: `kij = -0.054`.

They also proposed `kij = -0.0243 Cn + 0.1894` from methane through decane. This correlation is **not extrapolated** to the Middle or Heavy equivalent carbon numbers.

Baseline v1 therefore uses:

| pair | CPA initial kij | status |
|---|---:|---|
| H2O–Gasoline | 0.044 | n-hexane CPA proxy prior |
| H2O–Diesel | -0.054 | n-decane direct CPA anchor only |
| H2O–Middle | 0.0 | neutral prior; must be fitted |
| H2O–Heavy | 0.0 | diagnostic only; fitting forbidden |
| all oil–oil pairs | 0.0 | same-oil screening baseline |

These values are independent of the PR parameter table.

## Cross-association

The repository CPA interface can represent explicit cross association, as demonstrated by the validated water–DME benchmark. No water–oil-lump cross-association parameters are enabled in baseline v1.

For aromatic hydrocarbons, Oliveira et al. showed that CPA can require a solvation/cross-association parameter in addition to `kij`. Because Middle and Heavy are mixtures of saturates/aromatics/resins/asphaltenes rather than pure aromatics, importing benzene/toluene `beta_cross` would be scientifically unjustified.

## Calibration gate

CPA receives its own gate in `cpa_parameters/cpa_entry_gate.csv`.

A CPA binary is not accepted because the flash converges. It must reproduce:

1. phase count/topology;
2. both phase compositions;
3. density/volume data;
4. phase-boundary/critical-locus behavior;
5. hold-out states.

Association parameters are not allowed to compensate for an incorrect physical `a0/b/c1` set, and `kij` is not allowed to compensate for missing Heavy polarity.

## Current status

The CPA parameter system is now internally consistent and machine-readable, but it is **not calibrated for reservoir use**.

The immediate CPA-specific next step is a zero-flow water+n-hexane / water+dodecane / chemistry-bracket regression using these CPA parameters, followed by separate sensitivity for Middle aromatic solvation and a blocked Heavy branch until suitable heavy-oil data are available.
