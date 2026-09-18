# H02 run #86 mesh-normalized convergence audit

Audit date: 2026-09-19.

## Scope

This audit covers GitHub workflow `SCW H02 native 2D runner` run #86
(run `35371391415`, job `105685945849`, artifact `10559925922`).

All physical parameters remain frozen. The only numerical change relative to
the accepted run49/run64 package is the preregistered mesh-normalized SNES
absolute convergence gate:

```
RMS_ATOL  = 2.795084971874737e-8
LINF_ATOL = 1.0e-6
RTOL      = 1.0e-8
STOL      = 1.0e-100
```

The RMS threshold is exactly the old 20x8 `L2_ATOL=1e-6` divided by
`sqrt(1280)`. It was registered before the fine-grid rerun.

## Full-residual diagnosis that justified the change

The exact 60-s-aligned 60x20 B reproduction confirmed the old failure was a
global-L2 mesh-size effect, not a single bad thermodynamic row.

At the terminal plateau:

- global equation count: 9600;
- scaled L2 norm: approximately 1.00e-6 to 1.05e-6;
- scaled RMS norm: approximately 1.02e-8 to 1.07e-8;
- scaled Linf norm: approximately 1.12e-7 to 1.13e-7;
- the largest row was usually `O_W_BLOCK_OIL_HEAVY`;
- that largest row was only about 1.13e-7.

No individual equation was itself near 1e-6. The fixed global L2 absolute
tolerance therefore became artificially stricter as the number of equations
grew.

## 20x8 revalidation

A/B/C all completed beyond 2 actual PVI under the new RMS+Linf gate and passed
the existing strict global H2O/Heavy mass audit.

Therefore:

```
20X8_MESH_NORMALIZED_REVALIDATION = PASS
```

The accepted mechanism trend remains consistent with the previous strict
20x8 baseline.

## 60x20 B and C solver completion

Both fully compositional cases now complete beyond 2 actual PVI.

B:

- final time: 0.0701388889 day;
- accepted / rejected internal steps: 3168 / 92;
- minimum accepted dt: 7.41697e-7 day;
- maximum accepted dt: 2.3148148e-5 day.

C:

- final time: 0.0701388889 day;
- accepted / rejected internal steps: 3150 / 80;
- minimum accepted dt: 5.42535e-7 day;
- maximum accepted dt: 2.3148148e-5 day.

Thus the earlier 174-s artificial obstruction is removed for both B and C.

However, the workflow step was marked FAIL because the post-run numerical
audit correctly rejected the fine-grid global mass closure.

## Fine-grid global mass audit

Existing registered limit:

```
strict relative global mass error <= 1e-6
```

Observed maxima:

| mode | H2O | Heavy | result |
|---|---:|---:|---|
| B | 8.18e-7 | 3.13e-6 | FAIL (Heavy) |
| C | 7.59e-7 | 3.28e-6 | FAIL (Heavy) |

Phase-by-component producer accounting remains internally closed:

- instantaneous phase-component rate closure: 0 to numerical precision;
- cumulative phase-component closure: O(1e-17 kg).

The Heavy error is a reservoir-equation convergence accuracy issue, not a
producer-ledger accounting mismatch.

The Heavy balance error grows mainly during the middle/late displacement
period and reaches approximately:

- B: -1.80e-7 kg;
- C: -1.89e-7 kg.

The 20x8 mesh-normalized runs remain much tighter, so the RMS+Linf gate alone
is not sufficient to guarantee the pre-existing trajectory-level mass gate on
the fine mesh.

## Provisional grid sensitivity from completed trajectories

These values are useful diagnostically but are NOT yet certified grid-
convergence results because the 60x20 Heavy mass gate failed.

### 20x8 -> 60x20

| PVI | mode | RF difference | cumulative Heavy difference | DeltaP difference |
|---:|---|---:|---:|---:|
| 1 | B | -1.233 percentage points | 2.591% | 2.589% |
| 1 | C | -1.416 percentage points | 2.916% | 2.351% |
| 2 | B | -1.108 percentage points | 1.991% | 0.549% |
| 2 | C | -1.289 percentage points | 2.262% | 0.145% |

The preregistered grid gates were:

- RF <= 0.5 percentage point;
- cumulative Heavy <= 1%;
- pressure difference <= 5%;
- viscosity <= 1%.

Therefore even if the current fine-grid mass defect were ignored, the B/C
trajectories would already fail the RF and cumulative-Heavy grid screens.

Combined with the earlier A result (approximately 1% cumulative-Heavy
sensitivity), this makes the aligned 36x12 ensemble a necessary next grid
level once fine-grid numerical mass closure is restored.

## Why simply tightening RMS after seeing the result is not acceptable

The current fine-grid Heavy error could likely be reduced by lowering the RMS
threshold by an empirical factor. That would be a post-hoc fit to the observed
mass error and is not adopted.

The next numerical correction should directly constrain the quantity that must
remain mesh independent: the global component mass defect rate.

For the H02 2-PVI horizon:

- effective PV = 7.5e-5 m3;
- reservoir injection rate = 2.5e-8 m3/s;
- nominal 2-PVI time = 6000 s;
- the smaller relevant final conserved mass scale is about 0.055 kg;
- the existing trajectory mass target is 1e-6 relative.

A worst-case constant signed global component defect compatible with that gate
is therefore approximately

```
0.055 kg * 1e-6 / 6000 s = 9.2e-12 kg/s.
```

A conservative rounded gate of `8e-12 kg/s` for the absolute global signed
residual sum of each conserved component remains below the existing total
mass-error budget even through the 6060-s output overshoot.

This threshold is derived from the already registered mass gate, physical
inventory scale and fixed PVI horizon; it is not fitted to the observed
3.1e-6/3.3e-6 failures.

## Current status

```
FULL_RESIDUAL_EQUATION_DIAGNOSTIC        = PASS
GLOBAL_L2_MESH_SCALING_IDENTIFIED        = PASS
20X8_MESH_NORMALIZED_REVALIDATION        = PASS
60X20_B_SOLVER_COMPLETION                = PASS
60X20_C_SOLVER_COMPLETION                = PASS
60X20_B_GLOBAL_MASS_GATE                 = FAIL
60X20_C_GLOBAL_MASS_GATE                 = FAIL
CERTIFIED_20X8_TO_60X20_GRID_GATE        = BLOCKED
36X12_ALIGNED_GRID                       = REQUIRED_AFTER_MASS_FIX
REAL_HEAVY_PHYSICAL_VALIDATION           = NOT_VALIDATED
```

## Next small step

Add an opt-in, component-wise **global signed mass-residual convergence gate**
to the SNES convergence test:

```
abs(sum_cells R_component_raw) <= 8e-12 kg/s
```

for every conserved component, in addition to the existing RMS and Linf
conditions.

Then:

1. revalidate 20x8 A/B/C;
2. rerun only 60x20 B/C;
3. require the existing <=1e-6 trajectory mass audit;
4. only after that apply the registered grid metrics;
5. then run the 36x12 aligned-well A/B/C ensemble because the provisional
   20x8->60x20 RF/Heavy differences already exceed the grid thresholds.

No EOS, Heavy, kij, viscosity, rock or well parameter is to be changed.
