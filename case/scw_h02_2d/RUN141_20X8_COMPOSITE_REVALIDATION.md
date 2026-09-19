# H02 run #141 — 20x8 composite-callback revalidation

Audit date: 2026-09-19.

## Purpose

Revalidate the complete 20x8 A/B/C H02 trajectories after fixing SNES
convergence-callback ownership. The production Natural path now has one final
composite convergence callback that applies:

1. PETSc standard convergence/divergence semantics;
2. RMS residual gate;
3. residual infinity-norm gate;
4. component-wise global signed mass-residual gate;
5. stagnation detection only if the state has not been accepted.

No physical parameter or numerical tolerance was changed for this audit.

## Run provenance

- workflow: `SCW H02 native 2D runner`
- run: `35414796249` (#141)
- job: `105821189932`
- artifact: `10575896437`
- artifact SHA-256:
  `59295c79f46a1929092ea500286f6395e638f7e4934871732bb16ef2dd54a7c1`
- branch head:
  `d058ab8288125d55db86fa0ba30f89e2d693a28f`

All 60x20 and 36x12 long-run workflow steps were explicitly disabled for this
run. This audit certifies only the requested 20x8 revalidation.

## Fixed numerical gates

```
RMS_ATOL         = 2.795084971874737e-8
LINF_ATOL        = 1.0e-6
GLOBAL_MASS_ATOL = 8.0e-12 kg/s per conserved component
SNES_RTOL        = 1.0e-8
SNES_STOL        = 1.0e-100
SNES_MAX_IT      = 40
dt_max           = 2 s
```

The trajectory-level component mass requirement remains:

```
max relative component mass error <= 1e-6
```

## Complete trajectory result

A/B/C all reached 2.02 actual PVI.

| mode | accepted steps | rejected steps | max H2O mass error | max Heavy mass error |
|---|---:|---:|---:|---:|
| A | 3030 | 0 | 1.984e-7 | 3.358e-8 |
| B | 3044 | 8 | 7.098e-8 | 9.481e-9 |
| C | 3045 | 8 | 4.709e-8 | 1.153e-8 |

All registered numerical gates pass.

## Matched-PVI current result

| PVI | mode | RF_H | cumulative Heavy [kg] | DeltaP [MPa] | Heavy-carrier viscosity [Pa s] |
|---:|---|---:|---:|---:|---:|
| 1 | A | 0.397423117809 | 0.022897043731 | 0.003570723288 | 0.002000000000 |
| 1 | B | 0.476041577855 | 0.027426549532 | 0.003755925410 | 0.001999999876 |
| 1 | C | 0.485707136091 | 0.027983418772 | 0.003501183855 | 0.001732452134 |
| 2 | A | 0.483707830331 | 0.027868231232 | 0.002461032138 | 0.002000000000 |
| 2 | B | 0.556468634734 | 0.032060255414 | 0.002507342594 | 0.001999999725 |
| 2 | C | 0.569886098497 | 0.032833286072 | 0.002347604774 | 0.001732450380 |

## Comparison with run #49 and run #86

The comparison reuses the already registered 20x8 numerical-sensitivity
limits:

- RF_H absolute difference <= 0.001;
- cumulative Heavy relative difference <= 0.002;
- DeltaP relative difference <= 0.01;
- producer Heavy-carrier viscosity relative difference <= 0.005;
- B-A and C-B RF increment absolute difference <= 0.001.

Every gate passes against both references.

Largest observed differences versus run #49:

- RF_H: 3.435e-8 absolute (A, 2 PVI);
- cumulative Heavy: 7.101e-8 relative (A, 2 PVI);
- DeltaP: 1.238e-7 relative (A, 1 PVI);
- viscosity: effectively zero.

For B/C specifically, RF differences are only about 5e-9 to 7e-9 absolute and
cumulative-Heavy differences are about 1e-8 relative.

Run #86 gives essentially the same comparison; its stored B/C matched-PVI
values are identical to run #49 at recorded precision.

## Mechanism increments

Current run #141:

| PVI | B-A RF increment | C-B RF increment |
|---:|---:|---:|
| 1 | 0.0786184600461 | 0.0096655582360 |
| 2 | 0.0727608044026 | 0.0134174637630 |

Absolute differences from run #49:

| PVI | B-A difference | C-B difference |
|---:|---:|---:|
| 1 | 2.827e-8 | 1.375e-9 |
| 2 | 4.091e-8 | 1.221e-9 |

These are many orders of magnitude inside the 0.001 registered mechanism
consistency limit.

## Interpretation

```
20X8_COMPLETE_2PVI                    = PASS
20X8_TRAJECTORY_MASS_GATE             = PASS
RUN49_BASELINE_CONSISTENCY             = PASS
RUN86_BASELINE_CONSISTENCY             = PASS
B_MINUS_A_TREND_STABILITY              = PASS
C_MINUS_B_TREND_STABILITY              = PASS
PRODUCER_VISCOSITY_TREND_STABILITY     = PASS
```

Therefore the convergence-callback ownership repair does not materially change
the established 20x8 conditional mechanism result. It changes which nonlinear
states are legally accepted, not the converged physical branch on this mesh.

This remains
`CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION`.

## Next small step

Run only complete 60x20 B/C to 2 actual PVI with the same single composite
callback and unchanged gates. Do not rerun the already certified 60x20 A.
Require each B/C trajectory to pass the existing <=1e-6 component mass audit
before using any 20x8->60x20 grid comparison.
