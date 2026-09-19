# SCW-H02 A — true no-transfer control

This case is the mechanism-A counterfactual for SCW-H02. It uses the legacy independent-water formulation with one EOS hydrocarbon component, `OIL_HEAVY`. Water has its own conservation equation and never enters the Heavy EOS composition; Heavy never enters the independent water inventory.

Geometry, rock, wells, rate target, BHP cap, producer BHP, temperature, Heavy pure-component PR properties, relative permeability and actual-PVI ledger are matched to H02 B/C. This is a numerical no-transfer control, not a claim that real SCW and Heavy are immiscible. Do not emulate A by changing the B/C H2O-Heavy `kij`.

A shares the H02 managed solver defaults (`SNES_ATOL=1e-6`, `SNES_RTOL=1e-8`, `SNES_STOL=1e-100`, `SNES_MAX_IT=40`) so A/B/C use comparable numerical acceptance criteria. The fully compositional trace-phase appearance-hold correction affects B/C active-set behavior; it does not alter A's independent-water equations.

The accepted-step ledger integrates actual injector reservoir volume, and `Time::targetPVI=2.0` is active. A safety horizon remains mandatory; exhausting the horizon is not a target-PVI certificate.

Smoke:

```sh
make case CASE=scw_h02_2d_a -j2
make run CASE=scw_h02_2d_a NP=1 RESULT_DIR=$PWD/h02-A-smoke \
  RUN_ARGS='-numSteps 2 -dt 0.0000011574074074074074 -adaptive_dt true'
```

Complete verification:

```sh
make run CASE=scw_h02_2d_a NP=1 RESULT_DIR=$PWD/h02-A \
  RUN_ARGS='-target_pvi 2.0 -numSteps 150 \
  -dt 0.0006944444444444445 -adaptive_dt true \
  -ksp_type preonly -pc_type lu'
```

GitHub H02 native run #49 completed the 20x8 A verification trajectory beyond 2 actual PVI with maximum strict relative errors of approximately `4.57e-7` for H2O and `8.19e-8` for Heavy. Full 60x20 2-PVI grid convergence is still not established.

All A results remain a conditional numerical counterfactual, not experimental validation of immiscibility.
