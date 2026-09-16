# H2O-CO2-nC10 New-PR/New-SW/New-CPA SCW benchmark

This case reuses the production fully-compositional three-phase implementation
from `h2o_co2_nc10_2d_benchmark` at 653.15 K and 28 MPa. The 60 x 20 x 1 grid,
rock, wells, injection rate, and intended 0.1 PVI endpoint are unchanged.
Producer BHP is 27.5 MPa and injector maximum BHP is 30 MPa.

Select the thermodynamic model with `-eos pr`, `-eos sw`, or `-eos cpa`.

## SCW preflight status

The requested uniform oil/water state with saturations 0.80/0/0.20 is not an
admissible converged state of the current New-PR/New-SW/New-CPA implementation
at 653.15 K and 28 MPa. Temperature continuation and multi-seed restricted
flash scans found no two-phase root at the target condition. A common stable
overall composition, H2O/CO2/nC10 = 0.20/0/0.80, was also tested: PR and CPA
select a single oil-rich phase while SW selects a single water-role phase.
The flow-property preflight then correctly rejects the 0.80 nC10 water-slot
composition because the locked IAPWS-Garcia/McBride-Wright aqueous closure only
allows 1e-4 non-H2O/non-CO2 mole fraction.

Consequently this directory is a validation case, not a completed production
run. Do not loosen the aqueous-domain limit or present partial output as an
accepted simulation. See the validation report in
`RESULTS_SCW_653K_28MPa_NEW_MODELS_VALIDATION` before choosing a revised
initial state or a newly calibrated high-temperature phase model.
