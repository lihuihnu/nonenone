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

Stage 1 gates structural flash/stability/material closure only and reports
composition parity as an observed metric.  A composition-error threshold is not
invented until the remaining convention differences (site representation,
volume shifts, boundary feed, and phase-role mapping) are audited against the
paper.
