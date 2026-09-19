# Certified 60x20 A reference from H02 run #64

This directory contains the matched-PVI summary for the already completed
60x20 A no-transfer control. It exists so later aligned-grid studies can compare
36x12 A against the certified 60x20 A result without rerunning A.

Provenance:

- workflow: SCW H02 native 2D runner
- run: 35367347498 (#64)
- job: 105672991164
- artifact: 10557228665
- source artifact SHA-256:
  `4799bee6046ce6b68a2e76bd9b31683fa7326494ebf9996c8bfb2624192c25b8`
- source files:
  `convergence/dt2_60x20/A/producer_composition.csv`
  and `well_history.csv`
- A completed beyond 2 actual PVI and passed the strict global mass audit.
- Values are bracketed-linear interpolations at exactly 1 and 2 actual PVI,
  produced with the same formulas as `analyze_h02.py`.

This is a numerical conditional-mechanism reference, not physical validation.
