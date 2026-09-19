# H02 84x28 B/C short-grid final audit

Audit date: 2026-09-19.

## Scope

This audit combines the completed 84x28 B short trajectory from GitHub Actions
run `35426503074` with the completed 84x28 C-only trajectory from run
`35429797847`.

No physical parameter, convergence threshold, or grid-screening threshold was
changed.

All conclusions remain:

`CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION`.

## Provenance

### B

- source run: `35426503074`
- artifact: `10580795486`
- artifact SHA-256:
  `c0808223bcb50a4569d9ccebe0bcceddc8eae0bac6817f9c4516bf341c754c7c`
- head: `99bd1047988a7a77f666092556d1bbdea8c879cc`

The combined B/C workflow was cancelled by its 70-minute wall-clock limit
after B had completed. The B result itself is complete and valid.

### C

- source run: `35429797847`
- workflow conclusion: `success`
- artifact: `10580948700`
- artifact SHA-256:
  `cdddb8c6ee9f8dd3f4f756500742e5fedf9e7decd0473a9d7a01693145024807`
- head: `05047462e56b3162dfdb3f439b98baa518945129`

## Completion and numerical health

Both 84x28 modes reached approximately 1.22 actual PVI.

### B mass gate

- H2O max relative component mass error: `2.649025e-8`
- OIL_HEAVY max relative component mass error: `3.869239e-8`

### C mass gate

- final actual PVI: `1.219999999999`
- H2O max relative component mass error: `2.488264e-8`
- OIL_HEAVY max relative component mass error: `3.704883e-8`

Both are comfortably below the registered `1e-6` trajectory mass limit.

C completed with:

- 2102 accepted internal steps;
- 158 rejected/retried attempts;
- minimum accepted dt: about `0.1026 s`;
- maximum accepted dt: `2.0 s`;
- median accepted-step SNES iterations: `2`;
- maximum accepted-step SNES iterations: `14`.

There is no timestep-collapse signature.

## Registered grid limits

The existing H02 grid-screening limits are:

- Heavy RF absolute difference <= `0.005` = 0.5 percentage point;
- cumulative Heavy relative difference <= `1%`;
- DeltaP relative difference <= `5%`;
- producer Heavy-carrier viscosity relative difference <= `1%`.

The formal preregistered comparison locations were 1 and 2 PVI. The 0.5- and
1.2-PVI rows below are additional diagnostics only; no new formal gate is
invented there.

## 60x20 -> 84x28 B

| PVI | RF diff [pp] | cumulative Heavy diff | DeltaP diff | viscosity diff |
|---:|---:|---:|---:|---:|
| 0.5 | 0.248646 | 0.659386% | 2.027949% | ~4.5e-8% |
| 1.0 | **0.311771** | **0.672343%** | **0.147519%** | **~1.8e-7%** |
| 1.2 | 0.266594 | 0.549352% | 0.316395% | ~1.5e-7% |

At 1 PVI:

`60X20_TO_84X28_B_1PVI_GRID_GATE = PASS`.

## 60x20 -> 84x28 C

Matched values:

| PVI | RF_H 60x20 | RF_H 84x28 | cumulative Heavy 60x20 [kg] | cumulative Heavy 84x28 [kg] | DeltaP 60x20 [MPa] | DeltaP 84x28 [MPa] |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.38031866 | 0.37718613 | 0.02191158 | 0.02173111 | 0.00538357 | 0.00530055 |
| 1.0 | 0.47154576 | 0.46761922 | 0.02716757 | 0.02694131 | 0.00341894 | 0.00342947 |
| 1.2 | 0.49414334 | 0.49078456 | 0.02846948 | 0.02827595 | 0.00307538 | 0.00307028 |

Differences:

| PVI | RF diff [pp] | cumulative Heavy diff | DeltaP diff | viscosity diff |
|---:|---:|---:|---:|---:|
| 0.5 | 0.313253 | 0.823660% | 1.542145% | 0.000058% |
| 1.0 | **0.392654** | **0.832695%** | **0.308015%** | **0.000031%** |
| 1.2 | 0.335878 | 0.679718% | 0.165610% | 0.000028% |

At 1 PVI all four existing grid limits pass:

`60X20_TO_84X28_C_1PVI_GRID_GATE = PASS`.

## Combined B/C conclusion

The earlier 36x12 -> 60x20 aligned comparison failed marginally at 1 PVI:

- B cumulative Heavy: 1.0133%;
- C RF: 0.5513 pp;
- C cumulative Heavy: 1.1557%.

The next aligned refinement, 60x20 -> 84x28, now passes all four registered
integral metrics for both B and C at 1 PVI.

Therefore the earlier marginal failure is consistent with coarse-grid
discretization sensitivity rather than a nonconverged physical mechanism or
well-model defect.

The integral Heavy-recovery observables show a clear refinement trend:

- B 1-PVI RF shift:
  - 36x12 -> 60x20: about 0.475 pp;
  - 60x20 -> 84x28: about 0.312 pp.
- C 1-PVI RF shift:
  - 36x12 -> 60x20: about 0.551 pp;
  - 60x20 -> 84x28: about 0.393 pp.

The cumulative-Heavy differences similarly fall inside the 1% screen on the
60x20 -> 84x28 pair.

## Local-field caveat

Passing integral observables does not imply pointwise field convergence.

At C, 1 PVI, the high-Heavy tail continues to sharpen:

- area with `zHeavy > 0.5`:
  - 60x20: about 4.08%;
  - 84x28: about 5.31%;
- `zHeavy > 0.7`:
  - 60x20: about 2.75%;
  - 84x28: about 4.29%;
- `zHeavy > 0.9`:
  - 60x20: about 1.42%;
  - 84x28: about 2.85%.

This remains consistent with the previously identified front-sharpening /
coarse-grid numerical-diffusion mechanism.

Therefore the current claim is limited to convergence of registered integral
H02 observables, not pointwise convergence of the compositional front.

## Whether to extend 84x28 to 2 PVI

### For qualitative mechanism conclusions

No additional grid run is required merely to support the conditional ordering
`B > A` and `C > B`. That ordering has remained stable across the aligned
meshes.

### For formal full-range grid certification

Yes. The preregistered grid protocol requires matched comparisons at both
1 and 2 PVI.

The current 84x28 runs stop at 1.22 PVI, so they establish the previously
problematic 1-PVI region but cannot formally certify the 2-PVI endpoint.

The clean next step is therefore:

1. extend/run 84x28 B and C to 2 PVI with exactly the same frozen physics and
   composite convergence gates;
2. compare 60x20 -> 84x28 at exactly 2 PVI;
3. do not rerun A unless an unexpected B/C well-discretization regression
   appears.

No threshold should be changed.
