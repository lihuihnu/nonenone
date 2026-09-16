# five_component_eos_tuned_compare

## Purpose

This case compares three thermodynamic parameter packages after each package is
independently fitted to the same frozen public-data suite, during the same
five-component, three-phase reservoir displacement:

```text
H2O / CO2 / CH4 / C2H6 / nC4H10
```

All runs use the same `20 x 20 x 5` grid, rock, wells, initial `P-T-z`, TPFA
transport, fully compositional O/G/W flash, relative permeability and adaptive
time policy. Unlike `three_eos_3d_compare`, every backend receives its own
optimized pair parameters, while the public records, split protocol, optimizer,
reservoir and numerics remain common. Consequently, this is a comparison of
**best-fitted model packages**, not a single-factor comparison of EOS algebra.

## Thermodynamic packages and common-data fit

| token | package | intended strength |
|---|---|---|
| `pr` | PR78 alpha branch, fixed pure data, independently fitted constant BIPs | non-associating cubic baseline |
| `sw` | PR-SW water alpha, 1 mol/kg NaCl, fitted aqueous-correlation offsets and dry BIPs | salt-water/gas mutual solubility |
| `cpa` | SRK-CPA, fixed 4C-water/B2-CO2 association package, independently fitted BIPs | explicit association physics |

The three independent BIP matrices and exact pure/cross-association values are
recorded in `case_config.hpp`. The fitted pairs use 201 public VLE, solubility,
and gas-density records with a frozen 149/52 train/validation split. No
parameter is fitted to the generated reservoir results. Pairs absent from the
public suite retain their source-backed prior values.

The full protocol, fitted values, blind-validation errors, flow results,
limitations, and standard citations are recorded in
`COMMON_DATA_COMPARISON.md`. The executable and frozen data are under
`tools/example/five_component_common_fit/`.

Primary sources:

- Peng, D.-Y.; Robinson, D. B. *Ind. Eng. Chem. Fundam.* **1976**, 15,
  59-64. DOI `10.1021/i160057a011`.
- Søreide, I.; Whitson, C. H. *Fluid Phase Equilib.* **1992**, 77, 217-240.
  DOI `10.1016/0378-3812(92)85105-H`.
- Chabab, S. et al. *Int. J. Greenhouse Gas Control* **2019**, 91, 102825.
  DOI `10.1016/j.ijggc.2019.102825`.
- Qvistgaard, D. et al. *Fluid Phase Equilib.* **2023**, 570, 113796.
  DOI `10.1016/j.fluid.2023.113796`.
- Tsivintzelis, I.; Kontogeorgis, G. M. *J. Supercrit. Fluids* **2015**, 104,
  29-39. DOI `10.1016/j.supflu.2015.05.015`.

## Reservoir and wells

```text
grid/domain      20 x 20 x 5 / 1000 x 600 x 50 m
background       Kx/Ky/Kz = 80/50/5 mD, phi = 0.19
straight channel Kx/Ky/Kz = 250/160/18 mD, phi = 0.24
middle baffle    Kx/Ky/Kz = 2/1/0.05 mD, phi = 0.10
4 x 4 window     Kx/Ky/Kz = 150/100/10 mD, phi = 0.22
CO2_INJ          (i,j)=(1,4), k=0..1, 100000 reference m3/day
PROD             (i,j)=(18,15), k=3..4, BHP=52 bar
```

The common initial state is `P=60 bar`, `T=305 K`, and
`z=[0.30,0.10,0.15,0.15,0.30]`. The default schedule is 120 common targets of
0.25 day, for a total of 30 days. Internal steps are adaptive and transactional.
Outer boundaries are no-flow; mass enters or leaves only through the two wells.

## Build and run

```bash
make case CASE=five_component_eos_tuned_compare -j
make prepare CASE=five_component_eos_tuned_compare

make run CASE=five_component_eos_tuned_compare NP=2 EOS=pr
make run CASE=five_component_eos_tuned_compare NP=2 EOS=sw
make run CASE=five_component_eos_tuned_compare NP=2 EOS=cpa
```

By default, outputs are separated under `results/pr`, `results/sw`, and
`results/cpa`. A bounded pilot can override `NUM_STEPS`, `DT`, and `RESULT_DIR`.

"Tuned" means that each EOS receives its own optimum found under the common
bounded regression protocol. It does not mean that one universal parameter set
is mathematically optimal for every pressure, temperature or composition.

The reproducible 2026-08-25 pilot/full-run results, the single-donor and general
multi-donor CPA site-solver benchmarks, and engineering interpretation are
recorded in `RUN_REPORT.md`.

The full PR/SW/CPA recomputation after integrating the general CPA Newton
backend, its physical CSV equivalence audit, and the refreshed bilingual
figures are recorded in `RECOMPUTED_COMPARISON.md`.

The later density-corrected flow comparison adds IAPWS-anchored H2O volume
translation to PR/SW and a Garcia-based aqueous CO2 correction to SW. Its
preserved 30-day results are in
`results/translated_density_30day_20260826/{pr,sw,cpa}`;
the four-model bilingual figures and report are in
`../five_component_eos_tuned_legacy_compare/`.

## Clapeyron initial-flash cross-check

`clapeyron_flash_compare.jl` reproduces the common 6 MPa, 305 K initial flash
with Clapeyron.jl 0.6.27 using the exact PR78 and sCPA parameters in this case.
The seeded three-phase checks agree with the production results to `3.1e-5`
(PR78) and double-precision roundoff (sCPA); an optional unseeded full-TPD
search independently recovers the sCPA state. Clapeyron has no native
Soreide-Whitson model, so SW is not replaced by a non-equivalent surrogate.

See `CLAPEYRON_FLASH_COMPARISON.md` for conditions, numerical tables,
limitations, commands, and generated CSV paths.

## Result figures

Generate the bilingual publication-style figures from any completed result set:

```bash
python plot_results.py
# Use --force only when intentionally replacing an existing figure set.

python plot_results.py \
  --results-root results/recomputed_f405a18 \
  --figures-root figures/recomputed_f405a18 \
  --force
```

The Chinese and English PNG/PDF figures, the consolidated plotting table, and the
provenance manifest are written under `figures/`. See `RESULT_COMPARISON.md` for
the figure-by-figure engineering interpretation and quantitative comparison.

## Acceptance criteria

- all three initial flashes converge with O+G+W present;
- every common target is accepted without reaching the minimum time step;
- component inventories and well sources remain finite;
- the final component mass-balance report contains no failed conservation gate;
- PR/SW/CPA outputs are compared at identical physical times;
- solver effort is reported separately from physical response.
