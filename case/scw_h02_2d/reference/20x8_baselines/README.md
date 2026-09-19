# Certified 20x8 matched-PVI references

These two CSV files are retained verbatim from previously accepted H02
verification artifacts so the composite-convergence callback repair can be
checked against the established 20x8 mechanism trend.

- run49.csv
  - workflow run: 35365221250 (#49)
  - job: 105665913757
  - artifact: 10555863554
  - artifact SHA-256:
    2665dcbdea4f9922804e0381655a54d56ab2c851d855b0d554c2404398cf5ad1
  - source: full/analysis/h02_mechanism_at_1_2_pvi.csv
- run86.csv
  - workflow run: 35371391415 (#86)
  - job: 105685945849
  - artifact: 10559925922
  - source: mesh_norm_20x8/analysis/h02_mechanism_at_1_2_pvi.csv

Run #49 and run #86 B/C matched-PVI values are identical at stored precision;
A differs only at sub-numerical interpolation level.

Comparison thresholds reuse the already preregistered 20x8 numerical
sensitivity gates rather than introducing new post-hoc limits:

- RF_H absolute difference <= 0.001 (=0.1 percentage point)
- cumulative Heavy relative difference <= 0.002 (=0.2%)
- DeltaP relative difference <= 0.01 (=1%)
- producer Heavy-carrier viscosity relative difference <= 0.005 (=0.5%)
- B-A and C-B RF increment absolute difference <= 0.001

These are numerical consistency gates only. They do not validate real Heavy
thermodynamics.
