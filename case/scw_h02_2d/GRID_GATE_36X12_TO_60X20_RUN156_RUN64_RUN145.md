# H02 matched-PVI 36x12 -> 60x20 aligned grid gate

Audit date: 2026-09-19.

## Scope

This audit executes only the preregistered matched-actual-PVI aligned-grid gate.
No trajectory was rerun and no physical or numerical threshold was changed.

Certified inputs:

- 36x12 A/B/C: run #156, full 2 PVI composite-callback revalidation;
- 60x20 A: certified run #64 reference;
- 60x20 B/C: run #145, full 2 PVI composite-callback revalidation.

The 36x12 and 60x20 grids both place the intended injector and producer
completion centers exactly, so this comparison removes the completion-location
mismatch present in the earlier 20x8 -> 60x20 screen. Peaceman-WI and cell-size
effects still change with grid resolution.

## Provenance

36x12 run #156:

- workflow run: `35418851587`
- artifact: `10576886350`
- artifact SHA-256:
  `3a63e2bdf3b0bf42febab95046c8879f22f67d2f2d9b4de27e11fed7da33e5ae`
- head: `714319e766203f7c9fa59692f79994fb6f7d2028`

60x20 A run #64:

- workflow run: `35367347498`
- artifact: `10557228665`
- artifact SHA-256:
  `4799bee6046ce6b68a2e76bd9b31683fa7326494ebf9996c8bfb2624192c25b8`
- head: `21a695f53032a190f389b71413ca5d89c9b7d3bb`

60x20 B/C run #145:

- workflow run: `35415814058`
- artifact: `10575934743`
- artifact SHA-256:
  `168389652f1c0be4a0fdf39beee27805f5cb058a2aba8fe07172638ae4a40c98`
- head: `af0dbc1d40e31aba3f217b10c3d20ad0c9842ef7`

For 36x12 and 60x20 B/C, matched-PVI values are read from the generated
`analysis/h02_mechanism_at_1_2_pvi.csv` tables. The run #64 fine-grid A
matched-PVI values are reconstructed from its accepted-step
`producer_composition.csv` and `well_history.csv` using the same
`analyze_h02.py` bracketed-linear interpolation and A-mode viscosity rule.

## Preregistered thresholds

At matched actual PVI = 1 and 2:

```
abs(RF_H_60x20 - RF_H_36x12) <= 0.005
                               = 0.5 percentage point

relative cumulative Heavy difference <= 0.01
relative DeltaP difference           <= 0.05
relative producer Heavy viscosity    <= 0.01
```

The relative difference uses the registered symmetric scale

```
abs(test - base) / max(abs(test), abs(base))
```

Water-carried Heavy remains a trace-only diagnostic with no hard relative gate.

## Grid-gate result

```
36X12_TO_60X20_ALIGNED_GRID_GATE = FAIL
```

The failure is localized to the 1-PVI comparison. All A/B/C rows pass every
registered metric at 2 PVI.

| PVI | mode | RF difference [pp] | cumulative Heavy difference | DeltaP difference | viscosity difference | result |
|---:|---|---:|---:|---:|---:|---|
| 1 | A | 0.123844 | 0.314226% | 0.676885% | 0.00000000% | PASS |
| 1 | B | 0.474668 | 1.013263% | 0.379839% | 0.00000016% | **FAIL: cumulative Heavy** |
| 1 | C | 0.551336 | 1.155698% | 0.532718% | 0.00004812% | **FAIL: RF, cumulative Heavy** |
| 2 | A | 0.111303 | 0.231736% | 0.432443% | 0.00000000% | PASS |
| 2 | B | 0.391836 | 0.713328% | 0.718800% | 0.00000003% | PASS |
| 2 | C | 0.463246 | 0.824830% | 0.844708% | 0.00003306% | PASS |

The B 1-PVI cumulative-Heavy failure is marginal:
`1.013263%` versus the preregistered `1.0%` limit. It is nevertheless a
formal failure; the threshold is not relaxed after inspection.

The C 1-PVI RF difference is `0.551336` percentage point versus the
`0.5`-percentage-point limit, and cumulative Heavy differs by
`1.155698%`.

All pressure-drop and producer-viscosity metrics pass comfortably.

## Improvement relative to 20x8 -> 60x20

The aligned pair is materially less grid-sensitive than the previous
20x8 -> 60x20 comparison.

At 1 PVI:

- A RF difference decreases from about 0.454 pp to 0.124 pp;
- B RF difference decreases from about 1.233 pp to 0.475 pp;
- C RF difference decreases from about 1.416 pp to 0.551 pp.

This supports the earlier diagnosis that the 20x8 comparison mixed interior
resolution with completion-location / well-discretization effects. It does not
establish full grid independence because the aligned preregistered gate still
fails.

## Mechanism-increment sensitivity

No hard mechanism-increment threshold was preregistered, so these remain
diagnostic only.

| PVI | mechanism | 36x12 RF increment | 60x20 RF increment | change [pp] |
|---:|---|---:|---:|---:|
| 1 | B-A | 7.433077 pp | 7.082253 pp | -0.350824 |
| 1 | C-B | 0.860423 pp | 0.783755 pp | -0.076668 |
| 2 | B-A | 6.900849 pp | 6.620316 pp | -0.280534 |
| 2 | C-B | 1.231946 pp | 1.160536 pp | -0.071410 |

The qualitative ordering remains `B > A` and `C > B` on both aligned
meshes. The fine mesh reduces the magnitude of both mechanism increments.

## Interpretation boundary

This result is a numerical spatial/well-discretization sensitivity result for
the frozen H02 conditional mechanism experiment.

It does **not** validate the provisional OIL_HEAVY physical identity,
H2O-Heavy PR `kij=0.30`, or the B/C mixture-viscosity laws against real Heavy
experimental data.

Therefore all H02 conclusions remain:

`CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION`.

## Next small step

Do not rerun trajectories and do not retune physics.

Use the already saved 1-PVI 36x12 and 60x20 field snapshots to localize the
remaining B/C discrepancy into:

1. displacement-front spatial resolution;
2. completion-cell / Peaceman-WI resolution effects;
3. local phase-topology and Heavy-flux differences near the producer.

This diagnostic should use existing artifacts only and should not introduce a
new acceptance threshold after seeing the result.
