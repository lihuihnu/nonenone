# 60x20x1 supercritical-water displacement

This case is a controlled copy of the 60x20x1 H2O-CO2-nC10 benchmark. It retains the
300 m x 100 m x 5 m domain, 60x20x1 grid, rock properties, relative
permeabilities, well locations, injection rate, and 0.1-PVI comparison horizon.

The only physical control changed for the requested experiment is the uniform
isothermal state: 653.15 K and 28 MPa. Both exceed the pure-water critical point
(647.096 K and 22.064 MPa). Producer BHP is 27.5 MPa and injector maximum BHP is
30 MPa, so the prescribed well-pressure controls also remain above the critical
pressure.

Above the water critical point, the current fully-compositional O/G/W flash does
not admit the original 0.80/0/0.20 equilibrium initialization. This supported
variant therefore follows the repository's SCW convention: IAPWS supercritical
water is an independent conserved mobile phase and CO2-nC10 uses Peng-Robinson
phase behavior.

`RUN_PLAN.json` records units, boundaries, resource bounds, output requirements,
and numerical acceptance checks. The run is a sensitivity case; completion and
mass closure do not by themselves establish grid or timestep convergence.
