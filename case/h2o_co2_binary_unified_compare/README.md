# Unified H2O-CO2 binary comparison

This experiment uses one bulk initial composition for all models:
`z(H2O)/z(CO2)/z(nC10) = 0.80/0.20/0`.  It compares Traditional, New-PR,
New-SW, and New-CPA on the same 60 x 20 x 1 grid at 333.15 K/5.16 MPa
and 653.15 K/28 MPa.  The requested high-temperature condition is above
the pure-water critical temperature; the actual phase state is still determined
by each model's compositional flash.

The Traditional formulation has an independent pure-water phase, so its initial
gas/water saturations are calculated from the same bulk mole fractions and the
condition-specific phase molar densities.  The three new models flash the bulk
composition directly.  nC10 is intentionally absent; this avoids importing the
old oil-bearing initial state into a condition where the locked aqueous closure
does not support dissolved nC10.

Run definitions are in the four sibling wrapper cases.  Formal results are under
`outputs/h2o_co2_binary_unified_60x20x1/formal`.  The plotting script writes a
new, self-contained result directory with source data, validation tables, and
PNG/PDF/SVG figures.
