# H02 time-step and grid convergence protocol

Registered before convergence results are inspected.

## Scope

This protocol freezes the H02 physical/model package already accepted by the
Water-onset audit:

- 653.15 K;
- producer BHP 28 MPa;
- pure-water injector reservoir rate 2.5e-8 m3/s, 30 MPa maximum BHP;
- phi=0.25, k=1e-12 m2, geometry 0.30x0.10x0.010 m;
- provisional OIL_HEAVY identity and pure properties unchanged;
- H2O-Heavy PR kij=0.30 unchanged;
- B/C viscosity closures unchanged;
- kr=S^2 unchanged;
- trace-phase appearance hold unchanged;
- SNES_MAX_IT=40, SNES_ATOL=1e-6, SNES_RTOL=1e-8,
  SNES_STOL=1e-100;
- adaptive retry and accepted-step actual-PVI ledger unchanged.

No parameter may be retuned between the runs.

All cases stop only after actual accepted PVI reaches 2.0. Comparisons use
bracketed interpolation at exactly 1 and 2 PVI, never nominal time or first
point above target.

## Ensembles

1. 20x8x1, maximum internal dt = 2 s.
2. 20x8x1, maximum internal dt = 1 s.
3. 60x20x1, maximum internal dt = 2 s.

A/B/C are run for every ensemble. The convergence workflow uses the same MPI
rank count for all three ensembles so grid differences are not mixed with a
serial/parallel comparison.

## Independent numerical validity

Every ensemble must first pass the existing H02 audit:

- actual PVI coverage through 2;
- exact phase-component producer ledger closure;
- A has zero Water-carried Heavy;
- H2O strict global relative mass error <=1e-6;
- Heavy strict global relative mass error <=1e-6.

A failed ensemble is not used to claim convergence.

## Pre-registered time-step screening gates

For 20x8 dt_max=2 s versus 1 s, at PVI=1 and 2 for every A/B/C case:

- Heavy recovery fraction absolute difference <=0.001
  (=0.1 percentage point);
- injector-producer pressure-drop relative difference <=1%;
- cumulative Heavy production relative difference <=0.2%;
- producer Heavy-flux-weighted viscosity relative difference <=0.5%.

## Pre-registered grid screening gates

For 20x8 versus 60x20 at dt_max=2 s, at PVI=1 and 2 for every A/B/C case:

- Heavy recovery fraction absolute difference <=0.005
  (=0.5 percentage point);
- injector-producer pressure-drop relative difference <=5%;
- cumulative Heavy production relative difference <=1%;
- producer Heavy-flux-weighted viscosity relative difference <=1%.

Water-slot Heavy remains a trace-scale diagnostic (about 1e-8 of total Heavy
production in the accepted 20x8 study). Its absolute cumulative mass and share
are reported, but its relative change is not a hard convergence metric because
division by a near-zero trace would be misleading.

B-A and C-B mechanism increments in RF, pressure drop and viscosity are also
reported at matched PVI. No extra threshold is invented after seeing them.

## Interpretation

Passing these gates establishes only numerical screening convergence for this
conditional H02 mechanism experiment. It does not validate the provisional
Heavy pseudo-component, H2O-Heavy equilibrium, or mixture viscosity against
experimental data.
