# SCW-H02 A: true no-transfer baseline

This case is the mechanism-A control for SCW-H02. It uses the legacy
independent-water formulation with one EOS hydrocarbon component
(OIL_HEAVY). Water has its own conservation equation and never enters the
Heavy EOS composition; Heavy never enters the independent water inventory.

The geometry, rock, wells, rate target, BHP cap, producer BHP, temperature,
Heavy pure-component PR parameters, relative permeability and actual-PVI
ledger match H02 B/C. This is a numerical counterfactual baseline, not a
claim that real SCW and Heavy are immiscible.

Do not emulate A by changing the B/C H2O-Heavy kij.

Smoke:
```sh
make case CASE=scw_h02_2d_a -j2
make run CASE=scw_h02_2d_a NP=1 RESULT_DIR=$PWD/h02-A-smoke \
  RUN_ARGS='-numSteps 2 -dt 0.0000011574074074074074 -adaptive_dt true'
```

Full baseline uses Time::targetPVI=2.0 and the accepted-step actual reservoir
injection ledger. A safety horizon is still required; reaching the maximum
number of fixed output intervals without reaching target PVI is not a 2-PVI
success.
