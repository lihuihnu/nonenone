# H02 run #116 reciprocal-face conservation audit

Audit date: 2026-09-19.

## Purpose

This audit tested the hypothesis that the 60x20 B/C trajectory Heavy mass
drift was caused by evaluating each reciprocal TPFA connection independently,
so that an interior-face flux might not cancel exactly between its two cells.

No physical parameter or nonlinear tolerance was changed.

## Run

- workflow: SCW H02 native 2D runner
- run: 35406957069 (#116)
- artifact: 10573177553
- artifact SHA-256:
  `e8e484c05c6e49af6e0dc0abb8fcfc4cca79f0885360f1253095acf2ff8ce3ba`
- diagnostic case: 60x20 B
- target: 1.2 actual PVI
- RMS/Linf/global-signed-mass gates unchanged

For every accepted internal step the diagnostic recomputed all owned-cell
outward internal-face fluxes, MPI-summed the signed reciprocal contributions,
and integrated that global imbalance with the same backward-Euler accepted-step
time interval.

## Result at 1.2 PVI

| component | trajectory balance error [kg] | cumulative reciprocal-face imbalance [kg] | max instantaneous face imbalance [kg/s] |
|---|---:|---:|---:|
| H2O | +2.40446e-8 | +1.08446e-18 | 3.32254e-20 |
| OIL_HEAVY | -1.17487e-7 | +8.09084e-19 | 7.72024e-20 |

Thus the Heavy trajectory error exceeds the integrated face-pair imbalance by
approximately eleven orders of magnitude.

The diagnostic identity

```
trajectory_error + cumulative_internal_face_imbalance
```

is numerically indistinguishable from the trajectory error itself.

## Conclusion

```
RECIPROCAL_FACE_FLUX_NONCONSERVATION = REJECTED AS ROOT CAUSE
```

The current reciprocal face calculation is already conservative to floating
point for this H02 geometry. Refactoring to a single physical-face flux with
explicit +F/-F residual insertion may still be an architectural improvement,
but it cannot explain or repair the observed ~0.1 microgram-to-0.2 microgram
fine-grid Heavy trajectory drift.

The next root-cause branch is accepted-state consistency:

1. inspect the residual vector on which SNES certified convergence;
2. before commitTimeStep(), force a residual reevaluation using the exact final
   solution and independent phase-state vector;
3. compare cached vs fresh global signed component mass residuals;
4. if the fresh residual violates the registered 8e-12 kg/s gate, move the
   fresh-state certification into the convergence path;
5. if cached and fresh residuals agree, audit previous accumulation and accepted
   well-source identity per internal step.

All H02 interpretation remains conditional and not real-Heavy validation.
