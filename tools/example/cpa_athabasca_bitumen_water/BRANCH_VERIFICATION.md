# Branch reproducibility is not aggregate-AAD parity

## Audited baseline

Run 29, commit 417882075afb98d5367c1a93ae745be133f73a1e, artifact 10537793449.
Downloaded ZIP SHA256: 657273fd08474b40114eceed6684fbf273fa0d3b2715d6e5d2521acdc19c7c13.

The Amani-feed independent AAD is 0.7611628855469 MPa, but the continued-path AAD is 0.5312667403448 MPa. Matching the paper's aggregate 0.771 MPa is therefore NOT a reproduction certificate.

At 583.0 K the Amani-feed independent pressure is 11.50709368672 MPa while the continued pressure is 10.65566875532 MPa: a 0.85142493140 MPa path discrepancy. At 593.0 K the discrepancy is about 0.496129 MPa. These are not rounding effects.

At 573.1 K / 9.52 MPa, water mass fraction 0.559, unrestricted flash gives oil-phase xH2O=0.6468638008492; cold restricted OW gives 0.7202753438565 and fails missing-phase stability. Its material closure is nevertheless ~1.7e-18. Thus material closure plus a two-phase label does not certify equilibrium.

## New focused test

branch_consistency_probe.cpp directly includes the existing benchmark factory, without copying or tuning parameters. Five focused states compare unrestricted flash, cold restricted OW, and seeded pressure-continuation OW. Every active phase records stored/selected/liquid/vapor compressibility roots, molar density, composition and phase fraction. An independent same-root fugacity residual and dimensionless Gibbs value are recomputed.

Where unrestricted equilibrium is a stable OW state, the gate requires the other paths to recover the same certified branch (maximum composition difference <=1e-6). At genuine three-phase states it does not demand a metastable OW branch be globally stable. Invalid global certificates remain failures. These are numerical consistency criteria, not literature-fitting tolerances.

Possible mechanisms include different stationary branches or phase labels losing liquid/vapor root identity when Gibbs-root selection is enabled. The probe distinguishes these possibilities; neither is yet asserted as the sole root cause. Missing-phase density/root evidence and branch continuity must be checked before assigning a WLV-WL boundary.

OIL_HEAVY parameters, Athabasca literature parameters, and production defaults remain unchanged. Full Figure-7 reproduction and flow promotion remain BLOCKED until branch consistency and pointwise reference comparison are established.
