# Branch reproducibility is not aggregate-AAD parity

## Audited baseline

Run 29, commit 417882075afb98d5367c1a93ae745be133f73a1e, artifact 10537793449.
Downloaded ZIP SHA256: 657273fd08474b40114eceed6684fbf273fa0d3b2715d6e5d2521acdc19c7c13.

The Amani-feed independent AAD is 0.7611628855469 MPa, but the continued-path AAD is 0.5312667403448 MPa. Matching the paper's aggregate 0.771 MPa is therefore NOT a reproduction certificate.

At 583.0 K the Amani-feed independent pressure is 11.50709368672 MPa while the continued pressure is 10.65566875532 MPa: a 0.85142493140 MPa path discrepancy. At 593.0 K the discrepancy is about 0.496129 MPa. These are not rounding effects.

At 573.1 K / 9.52 MPa, water mass fraction 0.559, unrestricted flash gives oil-phase xH2O=0.6468638008492; cold restricted OW gives 0.7202753438565 and fails missing-phase stability. Its material closure is nevertheless ~1.7e-18. Thus material closure plus a two-phase label does not certify equilibrium.

Run 38, commit 8fa0724d7e8f40d141a5d876687b7fc96f76955e, artifact 10538775761, ZIP SHA256 8c575f243d58cba6d530246d1db1ed13faa2fe063b75f03fdb157d912e6f7d4a, additionally reports the 593.1 K / 12.77 MPa calibration audit: independent prediction 13.175295 MPa versus continued prediction 12.690592 MPa. This remains an audit, not a refit.

## Focused test and acceptance scope

branch_consistency_probe.cpp directly includes the existing benchmark factory, without copying or tuning parameters. Five focused states compare unrestricted flash, cold restricted OW, and seeded pressure-continuation OW. Every active phase records stored/selected/liquid/vapor compressibility roots, molar density, composition and phase fraction. An independent same-root fugacity residual and dimensionless Gibbs value are recomputed.

Where unrestricted equilibrium is a stable OW state, the gate requires the other paths to recover the same certified branch (maximum composition difference <=1e-6). When unrestricted flash has three phases, it does not demand a metastable OW branch be globally stable. Invalid global certificates remain failures. These are numerical consistency criteria, not literature-fitting tolerances. A three-phase slot assignment alone does not prove WLV topology.

## Confirmed result: branch-consistency run 1

Test commit: db093cafdfffbf31783ecb6986a780d5aa2aae03.
Workflow run: 35324062135; job: 105532919462.
Artifact: 10538631092; ZIP SHA256: 41d5cb20ef6f8a47c91e4a4c28519c0dd33efce40ec768685e66333d042a7706.
Compilation and artifact upload passed. The numerical branch-consistency gate failed: 3/5 focused states pass, 2/5 fail. This is a scientific/numerical failure, not an infrastructure failure.

This run preserves concurrent commit 994b863e2ee91db77a5fa64e1fcd4a815217be77, which changes the Athabasca solvation-beta interpretation. Consequently its numerical outputs must not be substituted into the older run-29/run-38 AAD calculations. The diagnostic commit itself changes no EOS factory or OIL_HEAVY parameter.

At 573.1 K / 9.52 MPa, water mass fraction 0.559:

| observable | cold restricted OW | continued OW / unrestricted |
|---|---:|---:|
| water fraction in oil | 0.720437 | 0.647174 |
| water fraction in Water slot | 0.984419 | 0.999949 |
| selected Water-slot Z | 0.665364 | 0.049812 |
| Water-slot molar density (mol/m3) | 3002.71 | 40108.61 |
| missing-phase stability | unstable | stable |
| dimensionless mixture Gibbs, G/(nRT), up to a common reference | -0.535014 | -0.572021 |

At the cold solution's own composition, Z_liquid=0.069967 while Z_vapor=0.665364, equal to its selected/stored Water-slot Z. Thus the slot named Water is demonstrably using the vapor-density root in this purported OW liquid branch. This is not an interpolation hypothesis: both roots, selected density, phase compositions, fugacity closure and stability were recalculated at the exact state.

Both cold and continued solutions satisfy material and active-phase fugacity closure to near roundoff. The cold solution nevertheless has higher total Gibbs and fails missing-phase stability. The same type of failure occurs for water mass fraction 0.441. The two failing cases have maximum full-composition path differences 0.07910029 and 0.07043332, respectively.

The 548.2 K control and 583 K / 11 MPa checks agree across paths. At 603.5 K / 15.32 MPa unrestricted flash finds three phases; OW-only solutions are unstable and are not silently promoted to globally stable equilibria. Physical identification of the third phase remains separate from its array slot name.

## Next correction, not yet applied

Separate the liquid-branch root constraint from unrestricted Gibbs-root selection. Preserve density-root identity while continuing the OW branch. Reject vapor-root aliases before interpreting a missing-phase stability flip as WLV-WL. Then rerun cold/warm/reference-state probes, verify all-component fugacity equality and total Gibbs ordering, and only afterwards regenerate Figure 7 pointwise.

The focused experiment establishes at least one concrete failure mechanism, not that every remaining phase-boundary discrepancy has the same cause. Do not tune kij, association parameters, or the desired AAD to conceal it.

Full Figure-7 reproduction and flow promotion remain BLOCKED. OIL_HEAVY parameters and production defaults are unchanged by this verification work.
