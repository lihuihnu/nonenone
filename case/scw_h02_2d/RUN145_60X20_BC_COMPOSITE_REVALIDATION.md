# H02 run #145 — 60x20 B/C composite-callback revalidation

Audit date: 2026-09-19.

## Scope

This audit revalidates only the full 60x20 B/C trajectories with the single
composite SNES convergence callback. The already certified 60x20 A reference
was not rebuilt or rerun. No 20x8 or 36x12 long-run grid study was executed.

## Provenance

- workflow: `SCW H02 native 2D runner`
- run: `35415814058` (#145)
- job: `105824146747`
- artifact: `10575934743`
- artifact SHA-256:
  `168389652f1c0be4a0fdf39beee27805f5cb058a2aba8fe07172638ae4a40c98`
- branch head:
  `af0dbc1d40e31aba3f217b10c3d20ad0c9842ef7`

## Fixed convergence gates

```
RMS_ATOL         = 2.795084971874737e-8
LINF_ATOL        = 1.0e-6
GLOBAL_MASS_ATOL = 8.0e-12 kg/s per conserved component
SNES_RTOL        = 1.0e-8
SNES_STOL        = 1.0e-100
SNES_MAX_IT      = 40
dt_max           = 2 s
```

Trajectory mass gate:

```
max relative component mass error <= 1e-6
```

## B acceptance

B completed 101 fixed output intervals, reaching 2.02 actual PVI.

- accepted / rejected internal steps: 3165 / 90
- nonlinear solves: 3255
- minimum accepted dt: 0.0640826 s
- maximum accepted dt: 2.0 s
- H2O max trajectory mass error: 2.787657e-8
- OIL_HEAVY max trajectory mass error: 5.071764e-8
- maximum Heavy absolute balance error: 2.922034e-9 kg

Result:

```
60X20_B_2PVI = PASS
60X20_B_MASS_GATE = PASS
```

## C acceptance

Because B passed, C was launched automatically and also completed 101 fixed
output intervals, reaching 2.02 actual PVI.

- accepted / rejected internal steps: 3150 / 80
- nonlinear solves: 3230
- minimum accepted dt: 0.046875 s
- maximum accepted dt: 2.0 s
- H2O max trajectory mass error: 3.587457e-8
- OIL_HEAVY max trajectory mass error: 3.071712e-8
- maximum Heavy absolute balance error: 1.769729e-9 kg

Result:

```
60X20_C_2PVI = PASS
60X20_C_MASS_GATE = PASS
```

## Improvement relative to run #86

Before the convergence-callback ownership fix, the same 60x20 trajectories
completed but failed the Heavy trajectory mass gate:

- B Heavy: ~3.13e-6 -> 5.07e-8, about 61.7x smaller
- C Heavy: ~3.28e-6 -> 3.07e-8, about 106.8x smaller

The fine-grid mass defect therefore tracks the callback-overwrite bug rather
than a change in physics or face-flux conservation.

## Matched-PVI values from run #145

| PVI | mode | RF_H | cumulative Heavy [kg] | DeltaP [MPa] | Heavy-carrier viscosity [Pa s] |
|---:|---|---:|---:|---:|---:|
| 1 | B | 0.463708 | 0.026716 | 0.003659 | 0.002000 |
| 1 | C | 0.471546 | 0.027168 | 0.003419 | 0.001732 |
| 2 | B | 0.545389 | 0.031422 | 0.002494 | 0.002000 |
| 2 | C | 0.556994 | 0.032091 | 0.002344 | 0.001732 |

Water-carried Heavy remains trace-only (~1e-8 to ~2e-8 cumulative fraction).

## Current certification status

```
20X8_COMPOSITE_REVALIDATION = PASS        # run141
60X20_A = PASS                            # certified run64 reference
60X20_B_COMPOSITE_REVALIDATION = PASS     # run145
60X20_C_COMPOSITE_REVALIDATION = PASS     # run145
60X20_BC_TRAJECTORY_MASS_GATE = PASS

20X8_TO_60X20_GRID_GATE = NOT YET EXECUTED AFTER FIX
36X12_ALIGNED_GRID = NOT YET EXECUTED
REAL_HEAVY_PHYSICAL_VALIDATION = NOT VALIDATED
```

## Next small step

Run only the preregistered matched-PVI 20x8 -> 60x20 grid comparison:

- A: run141 20x8 A versus certified run64 60x20 A reference;
- B/C: run141 20x8 B/C versus run145 60x20 B/C.

Do not change any grid-gate threshold. Record PASS/FAIL first; only after that
launch the already planned aligned 36x12 A/B/C study.
