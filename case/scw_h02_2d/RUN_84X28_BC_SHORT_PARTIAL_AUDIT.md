# H02 run #1 — 84x28 B/C short-grid partial audit

Audit date: 2026-09-19.

## Scope

This audit inspects GitHub Actions run `35426503074` for the dedicated
84x28 B/C short-grid study. The physical package and composite SNES gates were
frozen relative to the certified 36x12 and 60x20 H02 studies.

The workflow conclusion is `cancelled`, not `success`: the 70-minute job
limit was reached while mode C was running. The always-upload artifact was
preserved and is usable for the completed B trajectory and the partial C
trajectory.

Artifact:

- id: `10580795486`
- SHA-256:
  `c0808223bcb50a4569d9ccebe0bcceddc8eae0bac6817f9c4516bf341c754c7c`
- head: `99bd1047988a7a77f666092556d1bbdea8c879cc`

All conclusions remain:

`CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION`.

## Completion status

### B

B completed the requested 1.2-PVI target and wrote through approximately
1.22 actual PVI.

### C

C was numerically healthy but was interrupted by the workflow time limit at
approximately 0.28 actual PVI. This is not a nonlinear failure.

The final recorded C steps remain ordinary accepted steps with 2-3 SNES
iterations and `dt=2 s`; no timestep-collapse signature is present.

Partial C mass errors through 0.28 PVI:

- H2O max relative component mass error: `1.61e-8`;
- Heavy max relative component mass error: `7.90e-9`.

Therefore the C result is **incomplete due to wall-clock budget**, not failed
physics or failed nonlinear convergence.

## B mass certification

84x28 B full trajectory through 1.22 PVI:

- H2O max relative component mass error: `2.6490e-8`;
- Heavy max relative component mass error: `3.8692e-8`;
- Heavy max absolute mass error: about `2.23e-9 kg`.

Result:

`84X28_B_SHORT_MASS_GATE = PASS`.

## B: 36x12 -> 60x20 -> 84x28 matched-PVI trend

### 0.5 PVI

| grid | RF_H | cumulative Heavy [kg] | DeltaP [MPa] |
|---|---:|---:|---:|
| 36x12 | 0.3857765950 | 0.02222604365 | 0.00578105959 |
| 60x20 | 0.3770877968 | 0.02172544924 | 0.00577689973 |
| 84x28 | 0.3746013333 | 0.02158219471 | 0.00565974716 |

60x20 -> 84x28 diagnostic differences:

- RF_H: `0.248646 pp`;
- cumulative Heavy: `0.659386%`;
- DeltaP: `2.027949%`;
- Heavy-carrier viscosity: negligible.

0.5 PVI is diagnostic only; no new acceptance threshold is introduced there.

### 1.0 PVI

| grid | RF_H | cumulative Heavy [kg] | DeltaP [MPa] |
|---|---:|---:|---:|
| 36x12 | 0.4684548854 | 0.02698945159 | 0.00367263545 |
| 60x20 | 0.4637082065 | 0.02671597752 | 0.00365868535 |
| 84x28 | 0.4605904968 | 0.02653635451 | 0.00366409059 |

Using the already registered grid limits:

- RF absolute difference limit: 0.5 pp;
- cumulative Heavy relative difference: 1%;
- DeltaP relative difference: 5%;
- Heavy-carrier viscosity relative difference: 1%.

60x20 -> 84x28 gives:

- RF_H difference: **0.311771 pp — PASS**;
- cumulative Heavy difference: **0.672343% — PASS**;
- DeltaP difference: **0.147519% — PASS**;
- viscosity difference: `~1.8e-7%` — PASS.

Therefore:

`60X20_TO_84X28_B_1PVI_GRID_GATE = PASS`.

### 1.2 PVI diagnostic

60x20 -> 84x28:

- RF_H difference: `0.266594 pp`;
- cumulative Heavy: `0.549352%`;
- DeltaP: `0.316395%`;
- viscosity: negligible.

The B integral quantities continue to contract with refinement.

## B local-field behavior

The integral gate passes, but the local compositional front is still sharpening.

At 1 PVI:

- `zHeavy > 0.5` area:
  - 36x12: 0%;
  - 60x20: 4.00%;
  - 84x28: 4.34%.
- `zHeavy > 0.7` area:
  - 36x12: 0%;
  - 60x20: 2.00%;
  - 84x28: 3.10%.
- `zHeavy > 0.9` area:
  - 36x12: 0%;
  - 60x20: 0%;
  - 84x28: 2.17%.

Thus the earlier diagnosis remains supported: grid refinement primarily
sharpens the Heavy/water compositional front. Registered integral observables
are converging faster than pointwise/local front structure.

## Current decision

For mode B, the 84x28 short-grid result is sufficient to show that the earlier
36x12 -> 60x20 1-PVI marginal failure was a coarse-grid effect: the next
refinement 60x20 -> 84x28 passes all registered integral grid metrics.

Mode C is still required before the B/C aligned-grid study can be certified.

A separate C-only 84x28 workflow is used rather than rerunning B.
