# H02 run #156 — aligned 36x12 A/B/C composite revalidation

Audit date: 2026-09-19.

## Scope

This audit runs only the aligned 36x12 A/B/C ensemble to 2 actual PVI using
the frozen H02 physical package and the single composite SNES convergence
callback. No 20x8 or 60x20 long trajectory was rerun and no grid comparison
was executed in this workflow.

## Provenance

- workflow: `SCW H02 native 2D runner`
- run: `35418851587` (#156)
- job: `105832596885`
- artifact: `10576886350`
- artifact SHA-256:
  `3a63e2bdf3b0bf42febab95046c8879f22f67d2f2d9b4de27e11fed7da33e5ae`
- branch head:
  `714319e766203f7c9fa59692f79994fb6f7d2028`

## Frozen convergence gates

```
RMS_ATOL         = 2.795084971874737e-8
LINF_ATOL        = 1.0e-6
GLOBAL_MASS_ATOL = 8.0e-12 kg/s per conserved component
SNES_RTOL        = 1.0e-8
SNES_STOL        = 1.0e-100
SNES_MAX_IT      = 40
dt_max           = 2 s
```

Trajectory component-mass gate:

```
max relative component mass error <= 1e-6
```

## Completion and timestep behavior

All A/B/C trajectories completed 101 output intervals and reached about
2.02 actual PVI.

| mode | final actual PVI | accepted steps | rejected steps | nonlinear solves | min accepted dt [s] | max accepted dt [s] |
|---|---:|---:|---:|---:|---:|---:|
| A | 2.0199999919 | 3030 | 0 | 3030 | 2.0000 | 2.0000 |
| B | 2.0199999993 | 3060 | 19 | 3079 | 0.1976 | 2.0000 |
| C | 2.0199999997 | 3063 | 23 | 3086 | 0.2188 | 2.0000 |

No timestep-collapse behavior is present.

## Numerical mass audit

| mode | H2O max relative error | Heavy max relative error | Heavy max absolute error [kg] | result |
|---|---:|---:|---:|---|
| A | 1.986379e-7 | 3.172426e-8 | 1.827754e-9 | PASS |
| B | 2.268879e-8 | 1.593218e-8 | 9.179128e-10 | PASS |
| C | 1.847211e-8 | 1.763373e-8 | 1.015946e-9 | PASS |

Phase-component instantaneous closure is zero to reported precision. B/C
cumulative phase-component closure remains O(1e-17 kg).

Therefore:

```
36X12_A_2PVI            = PASS
36X12_B_2PVI            = PASS
36X12_C_2PVI            = PASS
36X12_TRAJECTORY_MASS   = PASS
```

## Matched-PVI values

| PVI | mode | RF_H | cumulative Heavy [kg] | DeltaP [MPa] | Heavy-carrier viscosity [Pa s] |
|---:|---|---:|---:|---:|---:|
| 1 | A | 0.394124113953 | 0.022706975684 | 0.003499652333 | 0.002000000000 |
| 1 | B | 0.468454885423 | 0.026989451586 | 0.003672635446 | 0.001999999876 |
| 1 | C | 0.477059120272 | 0.027485173986 | 0.003437220688 | 0.001732451295 |
| 2 | A | 0.480298648018 | 0.027671815389 | 0.002431549039 | 0.002000000000 |
| 2 | B | 0.549307141382 | 0.031647654790 | 0.002511636746 | 0.001999999738 |
| 2 | C | 0.561626603658 | 0.032357425444 | 0.002364230752 | 0.001732449932 |

Water-carried Heavy remains trace only:

- B: ~1.03e-8 at 1 PVI and ~2.29e-8 at 2 PVI;
- C: ~9.89e-9 at 1 PVI and ~2.21e-8 at 2 PVI.

## Interpretation boundary

This confirms numerical completeness and mass closure on the aligned 36x12
grid. It does not yet certify grid convergence. The 36x12->60x20 matched-PVI
comparison is intentionally deferred to a separate step.

All H02 conclusions remain:

`CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION`.

## Next small step

Execute only the preregistered matched-PVI 36x12 -> 60x20 grid gate using:

- 36x12 A/B/C: run #156;
- 60x20 A: certified run #64 reference;
- 60x20 B/C: run #145.

Do not rerun any trajectory and do not change any threshold.
