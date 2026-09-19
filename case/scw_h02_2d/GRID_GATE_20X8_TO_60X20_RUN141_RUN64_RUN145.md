# H02 matched-PVI 20x8 -> 60x20 grid gate

Audit date: 2026-09-19.

## Scope

This audit executes only the preregistered matched-actual-PVI grid gate. No new
trajectory was run.

Certified inputs:

- 20x8 A/B/C: run #141, composite callback, full 2 PVI revalidation;
- 60x20 A: certified run #64 reference;
- 60x20 B/C: run #145, composite callback, full 2 PVI revalidation.

All three source ensembles independently passed their applicable strict
trajectory mass audits before this comparison.

## Preregistered grid thresholds

At matched actual PVI = 1 and 2:

```
abs(RF_H_60x20 - RF_H_20x8) <= 0.005
                               = 0.5 percentage point

relative cumulative Heavy difference <= 0.01
relative DeltaP difference           <= 0.05
relative producer Heavy viscosity    <= 0.01
```

Water-carried Heavy is trace-reported only and has no hard relative gate.

No threshold was changed after observing the results.

## Grid-gate result

```
20X8_TO_60X20_GRID_GATE = FAIL
```

### Detailed matched-PVI differences

| PVI | mode | RF difference [pp] | cumulative Heavy difference | DeltaP difference | viscosity difference | result |
|---:|---|---:|---:|---:|---:|---|
| 1 | A | 0.4537 | 1.1417% | 2.6538% | 0 | FAIL: Heavy |
| 1 | B | 1.2333 | 2.5908% | 2.5890% | ~0 | FAIL: RF, Heavy |
| 1 | C | 1.4161 | 2.9156% | 2.3499% | 0.000097% | FAIL: RF, Heavy |
| 2 | A | 0.4522 | 0.9349% | 1.6253% | 0 | PASS |
| 2 | B | 1.1080 | 1.9911% | 0.5488% | ~0 | FAIL: RF, Heavy |
| 2 | C | 1.2892 | 2.2622% | 0.1425% | 0.000059% | FAIL: RF, Heavy |

All DeltaP and producer-viscosity gates pass.

A is only marginally outside the Heavy gate at 1 PVI, consistent with the
earlier run64 partial-grid audit. B/C show substantially larger spatial
sensitivity in RF and cumulative Heavy.

The fine grid consistently predicts lower RF and lower produced Heavy than the
20x8 grid:

- B RF shift: -1.233 pp at 1 PVI, -1.108 pp at 2 PVI;
- C RF shift: -1.416 pp at 1 PVI, -1.289 pp at 2 PVI.

## Mechanism-increment sensitivity

No hard mechanism-increment threshold was preregistered for the grid study, so
these values are diagnostic only.

| PVI | mechanism | 20x8 RF increment | 60x20 RF increment | change [pp] |
|---:|---|---:|---:|---:|
| 1 | B-A | 7.8618 pp | 7.0823 pp | -0.7796 |
| 1 | C-B | 0.9666 pp | 0.7838 pp | -0.1828 |
| 2 | B-A | 7.2761 pp | 6.6203 pp | -0.6558 |
| 2 | C-B | 1.3417 pp | 1.1605 pp | -0.1812 |

The qualitative ordering remains B > A and C > B on both meshes, but the
magnitude of both increments decreases on 60x20.

## Interpretation

This is a genuine spatial-discretization/well-discretization sensitivity result,
not a nonlinear-convergence artifact:

- run141 20x8 mass closure passed;
- run64 60x20 A mass closure passed;
- run145 60x20 B/C mass closure passed after the composite-callback fix.

The 20x8 mesh is therefore suitable as a low-cost screening grid but does not
satisfy the preregistered grid-independence screen for the B/C mechanism study.

The known well-discretization caveat remains important: 20x8 does not represent
the intended completion centers exactly, while 36x12 and 60x20 do. Therefore
20x8->60x20 mixes interior-grid refinement with completion-location and
Peaceman-WI changes.

## Next small step

Run the already planned aligned 36x12 A/B/C ensemble with the same frozen
physics and single composite convergence callback, then emphasize
36x12->60x20 as the cleaner spatial-convergence pair.

Do not change any physical parameter or grid-gate threshold.
