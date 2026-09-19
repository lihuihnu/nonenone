# H02 run #64 time/grid convergence audit

Audit date: 2026-09-19.

## Scope

This audit inspects GitHub workflow `SCW H02 native 2D runner` run #64
(run `35367347498`, job `105672991164`, artifact `10557228665`).

The physical/model package is frozen relative to the accepted Water-onset
audit. No EOS, Heavy, H2O-Heavy kij, viscosity, relative-permeability, rock or
well-control parameter is changed here.

## Run status

- trace Water appearance regression: PASS
- 60x20 two-step smoke: PASS
- strict first-60-s onset gate: PASS
- 20x8, dt_max=2 s, complete >2 PVI: PASS
- 20x8, dt_max=1 s, complete >2 PVI: PASS
- 60x20, dt_max=2 s:
  - A: complete >2 PVI
  - B: FAIL at about 174 s (~0.058 PVI)
  - C: not started because B failed
- final registered grid-convergence gate: not reached

The workflow conclusion is therefore FAIL, but the time-step convergence
evidence is complete and usable.

## Time-step convergence: 20x8, 2 s versus 1 s

All A/B/C runs independently pass the strict global H2O/Heavy mass gate.

Maximum strict relative mass errors for the 1-s ensemble:

| mode | H2O | Heavy |
|---|---:|---:|
| A | 1.78e-7 | 6.22e-8 |
| B | 1.21e-7 | 6.06e-8 |
| C | 1.17e-7 | 6.17e-8 |

At exactly matched actual PVI:

| PVI | mode | RF difference (percentage point) | DeltaP relative difference | cumulative Heavy relative difference | viscosity relative difference |
|---:|---|---:|---:|---:|---:|
| 1 | A | +0.00408 | 0.0229% | 0.0103% | 0 |
| 1 | B | -0.00473 | 0.00229% | 0.00995% | ~0 |
| 1 | C | -0.00633 | 0.000102% | 0.0130% | ~0 |
| 2 | A | +0.00316 | 0.0143% | 0.00654% | 0 |
| 2 | B | -0.000366 | 0.00178% | 0.000658% | ~0 |
| 2 | C | -0.00106 | 0.00701% | 0.00187% | ~0 |

These are all well inside the preregistered time-step gates:

- RF <= 0.1 percentage point;
- DeltaP <= 1%;
- cumulative Heavy <= 0.2%;
- viscosity <= 0.5%.

Therefore:

`20X8_TIME_STEP_CONVERGENCE = PASS`

The mechanism increments are also stable. At 2 PVI:

- B-A RF increment:
  - dt_max=2 s: +7.27608 percentage points
  - dt_max=1 s: +7.27255 percentage points
- C-B RF increment:
  - dt_max=2 s: +1.34175 percentage points
  - dt_max=1 s: +1.34105 percentage points

The dt=1 s runs do not reveal a materially different mechanism trend.

## 60x20 A result: useful partial grid evidence

A completed the full fine-grid trajectory with strict mass closure:

- maximum Heavy relative mass error: 4.69e-8
- maximum H2O relative mass error: 3.01e-7

Compared with 20x8 at matched PVI:

| PVI | RF 20x8 | RF 60x20 | RF abs difference | cumulative Heavy relative difference | DeltaP relative difference |
|---:|---:|---:|---:|---:|---:|
| 1 | 39.7423% | 39.2886% | 0.4537 percentage point | 1.1417% | 2.6538% |
| 2 | 48.3708% | 47.9186% | 0.4522 percentage point | 0.9349% | 1.6253% |

The preregistered grid limits were RF <=0.5 percentage point, cumulative Heavy
<=1%, DeltaP <=5%. A therefore already fails the cumulative-Heavy grid gate at
1 PVI (1.1417%), although the other A metrics pass.

Thus even before B/C are available:

`20X8_TO_60X20_GRID_CONVERGENCE = NOT_ESTABLISHED`

The result is close to the registered screening limits, not catastrophically
grid-sensitive.

## 60x20 B failure

B advances only to about 173.9 s (~0.058 PVI). The adaptive solver repeatedly
cuts the internal timestep and eventually reaches the configured minimum.

The final nonlinear plateau is:

- initial scaled SNES residual: ~1.12013e-1
- final scaled residual: ~1.001787e-6
- criterion: SNES_ATOL=1e-6
- reason: DIVERGED_LOCAL_MIN

The important diagnostic is that the conservation equations are already
closed at this failure. The failure diagnostic reports maximum local component
mass residual of order 1e-14 kg/s; the representative cell is an ordinary
Oil+Water cell with approximately:

- So = 0.3563
- Sg = 0
- Sw = 0.6437

and has no well source.

Therefore this is not the previous trace-Water appearance defect, and it is not
evidence of global mass leakage. The remaining ~1e-6 SNES norm must come from a
non-conservation equation block (phase-equilibrium/fugacity, volume closure or
well-control elsewhere). The current failure diagnostic does not identify that
row, so its exact identity is not yet certified.

Reducing timestep is not an effective remedy: adaptive retry already drives dt
from ordinary values to O(1e-5 s) and the residual plateau remains essentially
unchanged.

## Well discretization caveat

The 20x8 and 60x20 grids do not represent the wells identically.

20x8 cell centers for the mapped completions are approximately:

- injector: (0.0075, 0.01875) m
- producer: (0.2925, 0.08125) m

60x20 maps exactly to the intended locations:

- injector: (0.0125, 0.0125) m
- producer: (0.2875, 0.0875) m

The Peaceman WI also changes because cell dimensions change. Consequently
20x8->60x20 sensitivity contains spatial discretization, completion-location
and WI effects. It should not be labeled a pure interior-flow truncation error.

A 36x12 grid is a useful aligned intermediate: both target well centers are
represented exactly on 36x12 and 60x20, so a 36x12->60x20 comparison can
separate much of the coarse-well-location artifact from genuine field-grid
sensitivity.

## Recommended next step

Do not relax SNES_ATOL after observing the 60x20 failure, and do not retune
physics.

First add one failed-SNES diagnostic that records the largest **scaled full
residual** by:

- global cell / well row,
- equation identity (H2O mass, Heavy mass, O-W fugacity block, O-G block,
  volume closure, well control),
- residual value,
- phase presence and compositions.

Re-run only the short 60x20 B segment through ~200 s.

Then:

1. If the ~1e-6 plateau is a fugacity/closure row with physically acceptable
   per-row error and the issue is the mesh dependence of PETSc's global 2-norm
   absolute tolerance, preregister a mesh-independent convergence criterion
   (for example RMS/global-DOF normalization plus a fixed per-row infinity
   bound) and rerun both 20x8 baselines before the full 60x20 study.
2. If one specific thermodynamic row is genuinely stuck, fix that row/state
   transition rather than loosening the global tolerance.
3. Once B and C complete on 60x20, apply the already registered grid gates.
4. Because A is already marginal at 1 PVI, add a 36x12 aligned-grid ensemble if
   the 20x8->60x20 gate fails. Prefer 36x12->60x20 as the cleaner spatial
   convergence pair.

All conclusions remain conditional mechanism results, not real-Heavy physical
validation.
