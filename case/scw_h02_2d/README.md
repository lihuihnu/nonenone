# SCW-H02 native 2D case — staged implementation

This directory is executable with the repository case manager, not a design-only deck. Uses the production StructuredGridCore, NaturalStructuredGridRuntime, PR76 flash, AD finite-volume residual/Jacobian, Peaceman wells and accepted-step mass/producer ledgers. No CPA change.

Parameters: 60x20x1, 0.30x0.10x0.010 m, phi .25, k 1e-12 m2, 653.15 K, pure Heavy initial inventory, outlet 28 MPa, pure-water reservoir-volume injection 2.5e-8 m3/s with 30 MPa max BHP, kij .30. Components are H2O/OIL_HEAVY, not squalane. Planar horizontal grid with nz=1 has zero intercell gravity difference. Existing kr=S^2 retained.

Implemented modes B and C share the exact same PR. B viscosity uses pure-branch plateaus and a smooth log bridge; C uses log mixing of phase MASS fractions. The optional FluidSystem callback is unset by default so other cases remain unchanged. Mode A is explicitly rejected: the current full-composition runner does not implement the specified immiscible constraint. It is not faked with a different BIP.

```
make case CASE=scw_h02_2d -j2
make run CASE=scw_h02_2d NP=1 RESULT_DIR=$PWD/h02-C RUN_ARGS='-h02_mode C -numSteps 2 -dt 0.0000011574074074074074 -adaptive_dt true -snes_max_it 20 -snes_atol 1e-8 -snes_rtol 1e-8 -ksp_type preonly -pc_type lu'
```

-dt is in DAYS (0.0000011574074074074074 day = 0.1 seconds). -h02_nx/-h02_ny set mesh; well physical coordinates/radius fixed and WI recomputed. Never infer 2 PVI from number of requested steps: inspect actual producer_composition.csv PVI. Default is a TWO-STEP SMOKE run. 2 PVI automatic stopping, saturation-step bounds, target-PVI field interpolation, per-carrier component CSV and VTK conversion remain to be wired. Native CSV field snapshots and component ledgers are enabled now. Do not present partial native execution as the complete five-run study or an experimentally validated prediction. Schema and manifests retain these limitations.
