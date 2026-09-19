# Property validation

This directory contains **validation-only** thermophysical-property data for the SCW kerogen pseudo-component case.

- `density/`: aqueous density implementation checks and independent hydrocarbon/proxy density hold-outs.
- `viscosity/`: intrinsic transport checks using reference density, plus sources reserved for later coupled EOS+transport checks.

Rules:

1. validation rows are never used to fit PR/CPA BIPs;
2. validation rows are never used to fit volume translation or viscosity parameters;
3. source-reported expanded uncertainties are not reinterpreted as 1-sigma errors;
4. proxy success does not promote a pseudo-component to experimentally validated;
5. Heavy remains blocked even if squalane benchmark performance is good.

See `../08_DENSITY_VISCOSITY_VALIDATION.md`.
