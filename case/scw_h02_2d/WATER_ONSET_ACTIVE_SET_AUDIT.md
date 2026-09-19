# H02 Water-phase onset active-set audit

Audit date: 2026-09-18.

## Scope

This audit addresses one numerical failure only: trace Water-phase appearance in
the fully compositional H2O + OIL_HEAVY H02 B/C flow cases.

No EOS pure-component parameter, H2O-Heavy `kij=0.30`, viscosity closure,
relative permeability, rock property, geometry, well control or physical phase
criterion was fitted or changed to obtain the result.

The result remains a conditional mechanism experiment and does not validate the
provisional generated-Heavy pseudo-component against experiment.

## Original symptom

The earlier 20x8 H02 trajectory could complete when PETSc accepted
`CONVERGED_SNORM_RELATIVE`, but strict global mass closure then failed.

The first-60-s audit isolated accepted B/C steps with scaled final residuals of
order `1e-2`:

- B: the 20-22 s accepted step reached about `4.27e-2`;
- C: the 38-40 s accepted step reached about `1.87e-2`.

Mass error changed by discrete jumps at these events rather than drifting
continuously. Disabling step-norm convergence and requiring a scaled absolute
residual below `1e-6` removed the mass defect, but the 1-s B control then
exposed a separate active-set obstruction near 60 s rather than falsely
accepting it.

The obstructed cell was not a well cell. Its state was approximately

```
Oil saturation   ~ 0.999999744
Gas saturation   = 0
Water saturation ~ 2.6e-7
```

with local H2O/Heavy conservation residuals of order `1e-6 kg/s`.

## Active-set root cause

The production fully compositional active-set used:

- `phaseBoundaryProbeSaturation = 1e-4`;
- restricted flash + missing-phase stability after trace-phase removal;
- a wider reappearance margin only for phases recorded as suppressed.

The state machine had a gap at certified phase appearance:

1. Water was absent.
2. Missing-phase stability certified Water as unstable/required.
3. Global P-T-z flash legitimately reintroduced Water.
4. Near the phase boundary that flash could return a physically positive
   equilibrium saturation of order `1e-7`.
5. `assignFlashResult()` removed the suppression bit for every active phase.
6. On the next line-search post-check, the active Water phase satisfied
   `Sw <= 1e-4` and was deleted again solely because it was still in the
   numerical probe band.
7. The next missing-phase stability test reintroduced it.

The nonlinear iteration therefore alternated between a reduced state and a
certified trace-Water state. The problem was not a missing equilibrium solution
and was not a reason to alter `kij` or the phase-boundary probe itself.

## Correction: trace-phase appearance hold

The suppression state now also carries a short-lived appearance hold.

When a missing phase is reintroduced by a certified global flash:

- if `0 < S <= phaseBoundaryProbeSaturation`, its hold bit is retained;
- the phase remains active on subsequent post-checks;
- its physical saturation is left exactly at the flash/Newton value: no
  residual-saturation floor and no inflation to `1e-4` or `1e-6`;
- once `S > phaseBoundaryProbeSaturation`, the hold clears automatically;
- `S <= 0` is never protected and the ordinary restricted-flash/stability
  deletion path remains available.

Thus the fix changes active-set continuation memory, not thermodynamic
equations or equilibrium criteria.

Relevant implementation:

- `models/include/natural/state/three_phase_equilibrium.hpp`
- `models/include/natural/state/phase_state_data.hpp`
- focused regression in `test/src/unit/three_phase_flash_test.cpp`

## Focused regression

The unit regression explicitly constructs a positive Water phase at
`Sw=2.5e-7` with appearance hold and requires:

1. the phase is not immediately deleted;
2. `Sw` is not raised or otherwise changed by the hold;
3. the hold remains while `0 < Sw <= 1e-4`;
4. the hold clears after `Sw` grows to `2e-4`;
5. a negative Water saturation still enters ordinary phase-removal logic.

This regression passed in GitHub H02 native run #49.

## Real H02 validation

Workflow:

- name: `SCW H02 native 2D runner`
- run: `35365221250` (#49)
- job: `105665913757`
- artifact: `10555863554`
- artifact SHA-256:
  `2665dcbdea4f9922804e0381655a54d56ab2c851d855b0d554c2404398cf5ad1`

The run passed all of the following:

1. focused trace-phase active-set unit regression;
2. real 60x20 A/B/C two-step smoke;
3. first-60-s 20x8 differential matrix;
4. strict complete 20x8 A/B/C trajectories beyond 2 actual PVI;
5. phase-component producer closure;
6. strict global H2O and Heavy mass gates.

### First 60 s

The strict `dt_max=1 s` + `SNES_STOL=1e-100` +
`SNES_ATOL=1e-6` controls that previously exposed the B onset obstruction now
complete without rejected steps:

| case | Heavy relative mass error | H2O relative mass error | max final residual |
|---|---:|---:|---:|
| B strict dt=1 s | 1.44e-9 | 3.22e-8 | 6.39e-7 |
| C strict dt=1 s | 1.12e-9 | 2.50e-8 | 5.25e-7 |

The non-strict historical controls still reproduce the old behavior and are
kept as diagnostics: B/C can still accept high residuals if
`CONVERGED_SNORM_RELATIVE` is deliberately left enabled.

### Complete 20x8 strict 2-PVI trajectories

The registered strict full run uses:

```
SNES_ATOL = 1e-6
SNES_RTOL = 1e-8
SNES_STOL = 1e-100
SNES_MAX_IT = 40
maximum internal dt = 2 s
```

B and C each required eight rejected/cut internal attempts over the complete
trajectory. This is desired behavior: onset difficulty is handled by adaptive
retry instead of accepting a high-residual state.

Maximum strict relative global mass errors:

| mode | H2O | Heavy |
|---|---:|---:|
| A | 4.57e-7 | 8.19e-8 |
| B | 6.95e-8 | 1.61e-8 |
| C | 4.54e-8 | 1.70e-8 |

All are below the registered `1e-6` numerical gate.

## Effect on the conditional mechanism result

Cleaning the nonlinear/active-set path does not materially change the earlier
mechanism trend.

At exactly matched 2 actual PVI (bracketed interpolation, no extrapolation):

| mode | Heavy RF | injector-producer BHP difference | Heavy-flux-weighted producer viscosity |
|---|---:|---:|---:|
| A | 48.370786% | 2.4610 kPa | 2.0000 mPa s |
| B | 55.646863% | 2.5073 kPa | 2.0000 mPa s |
| C | 56.988609% | 2.3476 kPa | 1.73245 mPa s |

Water-slot Heavy remains a trace diagnostic:

- B: about `7.19e-10 kg` cumulative at 2 PVI,
  `2.24e-8` of total Heavy production;
- C: about `7.11e-10 kg`, `2.17e-8` of total Heavy production.

These values do not establish experimentally validated SCW extraction of Heavy.

## Production numerical policy

H02 managed launcher defaults are now registered in `case/Makefile` so the
mass-safe convergence policy is not CI-only:

```
SNES_MAX_IT = 40
SNES_ATOL   = 1e-6
SNES_RTOL   = 1e-8
SNES_STOL   = 1e-100
```

A uses the same acceptance tolerances for A/B/C comparability.

## Current gate state

```
TRACE_WATER_APPEARANCE_UNIT_REGRESSION = PASS
FIRST_60S_STRICT_ONSET_CROSSING         = PASS
20X8_STRICT_2PVI_TRAJECTORY             = PASS
GLOBAL_MASS_NUMERICAL_GATE              = PASS
PHASE_COMPONENT_LEDGER                  = PASS
FULL_60X20_2PVI_GRID_CONVERGENCE        = NOT_ESTABLISHED
REAL_HEAVY_PHYSICAL_VALIDATION           = NOT_VALIDATED
```

The next distinct task is time/grid convergence: compare strict 20x8
`dt_max=2 s` versus `1 s`, then complete the strict 60x20 trajectory at
matched actual PVI. That work must not alter the physical parameter set.
