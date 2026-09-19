# Athabasca bitumen + water CPA external proxy benchmark

This executable reproduces **Jia & Okuno (2018), Case 1** with the production
SRK-CPA backend.  It exists to test implementation capability for:

- standard 4C water;
- 4C self-associating asphaltene;
- water/asphaltene CR-1 cross association;
- non-self-associating pseudo-components that carry cross-association-only
  solvation sites;
- high-pressure bitumen/water Oil+Water(/Vapor) flash.

It does **not** calibrate or modify the kerogen `OIL_HEAVY` pseudo-component.
All Athabasca parameters are literature-benchmark-only.

Reference water-solubility rows are read from:

`case/scw_kerogen_lumped_flow_study/binary_pr_calibration/amani2013_athabasca_water_proxy.csv`

Stage 1 reproduces the **Oil+Water liquid branch** at all eight Jia & Okuno
Table 5 experimental WLV-WL transition states.  The hard gate requires a finite
O+W restricted solution with material closure <=1e-8 at every row.  This gate
passes 8/8 in the production CPA backend.

The unrestricted equilibrium result is deliberately reported separately: it
passes 7/8 rows.  At 603.5 K / 15.32 MPa the O+W branch exists and closes
material balance, but is unstable to an incipient missing phase.  That is not
used to fail Table 5 because those rows are experimental WLV-WL transition
points; a model phase-boundary offset can put the exact experimental T/P on the
neighboring topology.

Current observed branch-composition metrics are:

- O+W branch MAE versus experiment: about 0.01002 mole fraction;
- O+W branch MAE versus the published Jia CPA values: about 0.00491;
- maximum O+W material-closure error: about 6.7e-16.

No composition-error acceptance tolerance is invented.  The next distinct
benchmark is the Jia/Amani **Figure 7 WLV-WL phase boundary**, where unrestricted
phase topology and transition pressure can be gated against the reported
boundary data rather than inferred from Table 5 branch compositions.
