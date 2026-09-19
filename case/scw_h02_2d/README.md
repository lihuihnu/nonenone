# SCW-H02 native 2D case — conditional mechanism experiment

This directory is executable with the repository case manager. It uses the production `StructuredGridCore`, `NaturalStructuredGridRuntime`, PR76 flash, AD finite-volume residual/Jacobian, Peaceman wells, accepted-step mass ledgers and exact phase-by-component producer accounting. No CPA parameter or Heavy physical parameter is changed by this case.

## Fixed physical/model assumptions

The registered case uses a 60x20x1 topology over 0.30x0.10x0.010 m, porosity 0.25, permeability 1e-12 m2, 653.15 K, pure `OIL_HEAVY` initial inventory, producer BHP 28 MPa, pure-water reservoir-volume injection 2.5e-8 m3/s with 30 MPa maximum injector BHP, and H2O-Heavy PR `kij=0.30`. The Heavy component is the provisional >500 C pseudo-fraction, not squalane. The planar horizontal grid has zero intercell gravity difference and the current relative permeability closure remains `kr=S^2`.

Modes B and C use the exact same PR equilibrium model. B freezes the composition-viscosity feedback through the registered reference closure; C uses the current log-mixing rule on phase mass fractions. These are conditional mechanism assumptions, not measured high-temperature Heavy mixture properties.

Mode A is implemented separately in `../scw_h02_2d_a` as a true independent-water/no-transfer counterfactual. It is not emulated by changing B/C `kij`.

## Time/PVI semantics

`-dt` is the mandatory output interval in **days**, not the maximum nonlinear transport step. The default H02 output spacing is 60 s while `Time::maximumDtDays` caps accepted internal steps at 2 s. PVI is integrated from the **actual accepted injector reservoir rate**. `targetPVI=2.0` is implemented; reaching a fixed-step safety horizon without reaching target PVI is not a 2-PVI success.

The H02 managed launcher now defaults to:

- `SNES_ATOL=1e-6` on the scaled residual;
- `SNES_RTOL=1e-8`;
- `SNES_STOL=1e-100` so a small Newton step cannot certify a high-residual state;
- `SNES_MAX_IT=40`.

These numerical defaults were required because the old step-norm criterion could accept trace Water-phase onset states with scaled residuals of order 1e-2. They do not change EOS, `kij`, viscosity, rock or well physics.

## Water-phase onset active-set correction

The fully compositional active-set now preserves a short-lived **appearance hold** when stability legitimately reintroduces a phase with positive equilibrium saturation still inside the `1e-4` phase-boundary probe band. The physical trace saturation is not raised to a floor. The hold clears once the phase grows above the probe band; a zero or negative saturation can still trigger ordinary restricted-flash/stability removal.

This removes the deterministic sequence

`missing Water -> stability says reappear -> flash gives Sw~1e-7 -> next post-check deletes Water`

that previously trapped B near `Sw~2.6e-7`.

See `WATER_ONSET_ACTIVE_SET_AUDIT.md`.

## Verified execution state

GitHub H02 native run #49 (workflow run `35365221250`, artifact `10555863554`) passed:

- focused trace-phase active-set unit regression;
- real 60x20 A/B/C two-step smoke;
- 20x8 first-60-s differential diagnosis, including 1 s internal-step controls;
- complete strict 20x8 A/B/C trajectories beyond 2 actual PVI;
- phase-component producer closure and global H2O/Heavy mass gates.

For the complete 20x8 verification trajectories, maximum strict relative mass errors were:

| mode | H2O | Heavy |
|---|---:|---:|
| A | 4.57e-7 | 8.19e-8 |
| B | 6.95e-8 | 1.61e-8 |
| C | 4.54e-8 | 1.70e-8 |

All are below the registered `1e-6` numerical gate. B/C required eight rejected/cut internal attempts over the full trajectory; the onset is therefore crossed through normal adaptive retry rather than by accepting a high-residual state.

The full 60x20 **2-PVI** grid study has not yet been run. The 60x20 result currently certified by CI is the short smoke only. Grid convergence is therefore not established.

## Example

```sh
make case CASE=scw_h02_2d -j2
make run CASE=scw_h02_2d NP=1 RESULT_DIR=$PWD/h02-C \
  RUN_ARGS='-h02_mode C -target_pvi 2.0 -numSteps 150 \
  -dt 0.0006944444444444445 -adaptive_dt true \
  -ksp_type preonly -pc_type lu'
```

Use `-h02_nx/-h02_ny` only for explicit verification-grid studies; well physical coordinates/radius remain fixed and the well index is recomputed. Always compare cases on actual PVI from `producer_composition.csv`.

Every H02 result remains:

`CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION`

A numerically closed trajectory does not promote the provisional Heavy equilibrium or mixture-viscosity assumptions to experimental validation.
