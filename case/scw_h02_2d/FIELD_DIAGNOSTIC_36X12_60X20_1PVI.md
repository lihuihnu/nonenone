# H02 1-PVI field diagnostic: 36x12 versus 60x20

Audit date: 2026-09-19.

## Scope

This diagnostic uses only already-saved H02 artifacts. No trajectory was rerun,
no physical parameter was changed, and no new acceptance threshold was
introduced.

Inputs:

- aligned 36x12 A/B/C: run #156;
- 60x20 A: certified run #64;
- 60x20 B/C: run #145.

The comparison focuses on the saved step-50 fields, which are at essentially
exactly 1 actual PVI. Existing step-25 fields at 0.5 PVI are used only to
identify when the cumulative Heavy gap was created.

All conclusions remain:

`CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION`.

## Question

The aligned 36x12 -> 60x20 grid gate failed only near 1 PVI:

- B cumulative Heavy: 1.013263% difference versus 1% limit;
- C RF_H: 0.551336 percentage point versus 0.5 pp limit;
- C cumulative Heavy: 1.155698% versus 1% limit.

At 2 PVI every A/B/C metric passed.

This diagnostic separates four possible causes:

1. Heavy/water-front spatial resolution;
2. phase-topology resolution;
3. producer-neighborhood Heavy flux;
4. completion-cell / Peaceman-WI discretization.

## Exact 1-PVI field alignment

Step 50 corresponds to 3000 s and approximately 1 actual PVI for every run:

| case | PVI |
|---|---:|
| 36x12 B | 0.999999999999467 |
| 60x20 B | 0.999999999999481 |
| 36x12 C | 1.000000000003261 |
| 60x20 C | 0.999999999999156 |

Therefore the field comparison does not mix PVI mismatch into the spatial
diagnostic.

The registered producer completion center is represented exactly on both
aligned grids:

- 36x12 producer input cell: 394;
- 60x20 producer input cell: 1077.

## 1. Front resolution: primary cause

A conservative common-area comparison was made on a 180x60 overlay. This is an
exact common subdivision of both uniform grids:

- every 36x12 cell maps to 5x5 common subcells;
- every 60x20 cell maps to 3x3 common subcells.

No interpolation or smoothing is required for the area comparison.

### Water-saturation sharpness

At 1 PVI:

| metric | B 36x12 | B 60x20 | C 36x12 | C 60x20 |
|---|---:|---:|---:|---:|
| mean Sw | 0.376464 | 0.374975 | 0.387824 | 0.385166 |
| Sw standard deviation | 0.153965 | 0.160940 | 0.160535 | 0.169286 |
| total-variation diagnostic | 0.197528 | 0.227189 | 0.204839 | 0.232560 |

The mean water content changes little, while the fine grid has substantially
larger variance and total variation. This is the signature expected when the
fine grid resolves a sharper displacement/compositional front and the coarser
grid numerically smears it.

The same tendency is already present at 0.5 PVI:

| metric | B 36x12 | B 60x20 | C 36x12 | C 60x20 |
|---|---:|---:|---:|---:|
| mean Sw | 0.287781 | 0.283224 | 0.294850 | 0.289755 |
| total-variation diagnostic | 0.198611 | 0.217857 | 0.204161 | 0.223622 |

### Heavy-composition tail

The sharper fine-grid front is even clearer in total Heavy mole fraction
(`zHeavy = 1 - zH2O`).

At 1 PVI:

- B:
  - 36x12 has **0%** of area with `zHeavy > 0.5`;
  - 60x20 retains **4.0%**;
  - 36x12 has **0%** with `zHeavy > 0.7`;
  - 60x20 retains **2.0%**.
- C:
  - `zHeavy > 0.5`: 2.083% versus 4.083%;
  - `zHeavy > 0.7`: 0% versus 2.75%;
  - `zHeavy > 0.9`: 0% versus 1.417%.

At 0.5 PVI the same pattern is stronger:

- B `zHeavy > 0.7`: 4.630% on 36x12 versus 7.583% on 60x20;
- C `zHeavy > 0.7`: 5.324% versus 9.083%.

Thus the 36x12 mesh removes the high-Heavy tail more rapidly through spatial
smearing. This directly explains why it transports/produces slightly more
Heavy before 1 PVI.

## 2. Phase topology: secondary manifestation of the same front-resolution issue

At 1 PVI the B/C phase masks contain Oil-only and Oil+Water states.

| metric | B 36x12 | B 60x20 | C 36x12 | C 60x20 |
|---|---:|---:|---:|---:|
| Oil-only area | 2.315% | 4.167% | 3.472% | 4.750% |
| common-area phase-mask mismatch | colspan | 3.704% | colspan | 3.278% |

The 60x20 grid preserves more Oil-only area, consistent with a sharper water
onset front. The topology mismatch is not concentrated at the producer:
less than about 3% of the mismatch area lies within 3 cm of the producer.

Most mask disagreement occurs in the upper boundary strip, where the 36x12
and 60x20 vertical cell sizes resolve the steep phase-onset layer differently.
This is a spatial-resolution effect rather than evidence of a different
thermodynamic branch.

## 3. Producer completion cell: nearly converged

At exactly 1 PVI:

### B

- producer-cell Sw:
  - 36x12: 0.2848767
  - 60x20: 0.2821714
  - absolute difference: 0.0027052
- oil-phase H2O mole fraction:
  - 36x12: 0.59747344
  - 60x20: 0.59747433

### C

- producer-cell Sw:
  - 36x12: 0.2937022
  - 60x20: 0.2914187
  - absolute difference: 0.0022835
- oil-phase H2O mole fraction:
  - 36x12: 0.59747369
  - 60x20: 0.59747453

The local equilibrium state at the completion cell is therefore almost the
same. The grid failure is not caused by a large local phase-state discrepancy
at the producer cell.

## 4. Peaceman WI: not the cause of the lower fine-grid cumulative RF

The producer total well index changes, as expected, when cell size changes:

- 36x12: `2.9775840735e-14`;
- 60x20: `3.9286198222e-14`;
- fine/coarse ratio: **1.319398** (+31.94%).

However the cell-to-BHP pressure drawdown compensates:

### B at 1 PVI

- completion-cell drawdown:
  - 36x12: 448.390 Pa
  - 60x20: 345.599 Pa
  - fine/coarse ratio: 0.77076
- `WI * drawdown` fine/coarse ratio: 1.01693
- total reservoir production rate fine/coarse: 1.00138

### C at 1 PVI

- drawdown:
  - 36x12: 419.178 Pa
  - 60x20: 321.640 Pa
  - fine/coarse ratio: 0.76731
- `WI * drawdown` fine/coarse ratio: 1.01239
- total reservoir production rate fine/coarse: 0.99989

So the approximately 32% WI increase mainly changes the local pressure
drawdown required by the fixed-BHP well. It does not create a comparable total
production-rate change.

The independent-water A control reaches the same conclusion. A experiences
the same grid/WI change, but at 1 PVI its aligned-grid differences are only:

- RF_H: 0.123844 pp;
- cumulative Heavy: 0.314226%.

Both pass comfortably. The much larger B/C sensitivity therefore cannot be
attributed primarily to WI.

## 5. Instantaneous Heavy flux at 1 PVI acts opposite to the cumulative gap

At 1 PVI the fine grid is **not** producing less Heavy instantaneously.

### B

- Heavy mass rate:
  - 36x12: `2.2335359e-6 kg/s`
  - 60x20: `2.2885725e-6 kg/s`
  - fine grid: **+2.464%**
- oil reservoir-volume rate: also +2.464%;
- total reservoir-volume rate: only +0.138%.

### C

- Heavy mass rate:
  - 36x12: `2.3513550e-6 kg/s`
  - 60x20: `2.3959072e-6 kg/s`
  - fine grid: **+1.895%**
- total reservoir-volume rate differs by only -0.011%.

Therefore the lower 60x20 cumulative RF at 1 PVI is a history effect. It was
created before 1 PVI and is already being reduced by 1 PVI.

## 6. When the cumulative difference is created

The largest coarse-minus-fine cumulative Heavy gap occurs at about **0.5 PVI**,
not at 1 PVI.

### B

- maximum gap at 0.5 PVI: `5.00594e-4 kg`;
- relative to the larger cumulative value: about **2.25%**;
- by 1 PVI the gap has already shrunk to `2.73474e-4 kg` (**1.013%**);
- by 2 PVI: `2.25752e-4 kg` (**0.713%**).

At 0.24 PVI the fine grid Heavy rate is about **4.65% lower** than the coarse
grid. By 0.74 PVI it is about **5.12% higher**, so the fine grid is already
catching up well before 1 PVI.

### C

- maximum gap at 0.5 PVI: `5.62968e-4 kg`;
- relative gap: about **2.50%**;
- 1 PVI: `3.17646e-4 kg` (**1.156%**);
- 2 PVI: `2.66894e-4 kg` (**0.825%**).

At 0.24 PVI the fine-grid Heavy rate is about **8.83% lower**; by 0.74 PVI it
is about **5.50% higher**.

This temporal behavior is exactly what is expected from a sharper fine-grid
front: less early numerical spreading/Heavy displacement, followed by a later,
sharper arrival and catch-up.

## Attribution

The evidence supports the following ordering.

### Primary

**Displacement/compositional-front spatial resolution.**

The 36x12 mesh is more diffusive. It smooths water and Heavy composition over a
larger region, removes the high-Heavy tail earlier, and produces too much Heavy
during the early/mid displacement period. The fine grid preserves a sharper
front and catches up after approximately 0.5 PVI.

### Secondary

**Phase-onset/topology resolution, especially near the upper boundary layer.**

This is consistent with the same front-resolution mechanism; it is not a
separate EOS-branch failure.

### Small / not causal for the observed sign

**Completion-cell / Peaceman-WI discretization.**

WI changes substantially, but pressure drawdown compensates and total well rate
barely changes. At 1 PVI the fine-grid instantaneous Heavy rate is actually
higher, opposite to the cumulative-RF deficit.

## Decision on another aligned grid

Two different goals must be separated.

### If the goal is qualitative H02 mechanism isolation

No additional grid is required to support the qualitative ordering:

- B > A in Heavy recovery;
- C > B in Heavy recovery;
- the ordering survives both 36x12 and 60x20;
- 2-PVI aligned-grid gates already pass.

The current H02 result is still conditional and not real-Heavy validation.

### If the goal is formal numerical grid-independence certification

One finer aligned grid is required because the preregistered 36x12 -> 60x20
gate formally fails at 1 PVI.

The next uniform aligned topology should be **84x28**.

For the fixed 0.30 m x 0.10 m geometry and well centers
(0.0125, 0.0125) m and (0.2875, 0.0875) m, the uniform aligned sequence is:

- 36x12;
- 60x20;
- **84x28**;
- 108x36; ...

84x28 therefore preserves exact completion-center alignment while increasing
cell count from 1200 to 2352, a modest factor of 1.96.

### Recommended staged execution

Do **not** rerun A first.

1. Run B/C only on 84x28 through about 1.1-1.2 actual PVI.
2. Compare 60x20 -> 84x28 at exactly 0.5 and 1.0 PVI using the already
   registered RF / cumulative-Heavy / DeltaP / viscosity metrics.
3. Inspect whether:
   - the early 0.5-PVI cumulative gap contracts;
   - 1-PVI RF and cumulative Heavy enter the existing thresholds;
   - the high-`zHeavy` tail and Oil-only area approach a stable limit.
4. Only if 60x20 -> 84x28 is satisfactory should the 84x28 run be extended to
   2 PVI for formal full-range certification.
5. A need not be repeated unless the B/C result suggests an unexpected
   well-discretization regression.

No threshold should be modified after seeing the 84x28 result.
