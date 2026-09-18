# Identified vapor-onset continuation versus oil composition spinodal

## Scope and execution provenance

This is a diagnostic, not a fitted phase-boundary curve and not a new production flash solver. Heavy parameters, the frozen Athabasca/Jia factory, production EOS equations, existing flash tolerances, production defaults and flow/reproduction gates are unchanged.

The local execution reused the exact production source from artifact 10541515337 (source checkout 303bda33d2ef5e67e6a10c7300f25ae6250e626b). Comparing its source revision 39e0af4fafce7b7dd8ca37acf683bb210b453947 with branch revision da01f2ed06face8e942e956f9fc383842147c257 showed only the preceding identity-audit README/analyze.py additions. The frozen factory Git blob was independently checked as f23ad37e7673e756274e36762ce7b3f43aa1c3ad.

New source commits: b782bf7aa0f1b0f39f5138d3fbd1b3f1e8909765 (main.cpp) and ea9ce3c5f16cbf589fe4d90f8229aefad4137533 (refine.cpp). Their remote Git blobs match the locally executed files exactly:

- main.cpp: 35c4c1cc28e78aaea738c7db051eb60e7659826d; SHA256 90561adced7b50a7043b0c06523d3c35348a91605473bf26551813ec63d29180.
- refine.cpp: d1fb9c0c325b43c56a08957159f49f24c6391cc3; SHA256 51ca183b030f466fb86ee7dc36c684544366ed34b399c924b5a08fa82f88c02e.

Both compiled locally with C++17, O2, Wall/Wextra/Wpedantic/Werror and completed. The separate read-only CI records its own checkout and outputs; local results below must not be relabelled as remote CI results.

## Equations and branch identity

Feed is the authoritative 44.1 wt% water / 55.9 wt% bitumen convention. Start from the nontrivial 583 K candidate identified by the production fixed-root OW flash/stability routines, then solve the incipient three-phase equations directly:

- mu_i(O) = mu_i(W) = mu_i(V), all five components;
- z_i = beta_O*x_i(O) + (1-beta_O)*x_i(W), with beta_V=0;
- all three compositions normalized; pressure is unknown at each prescribed T.

Three sets of four log composition ratios, logit(beta_O), and log(P/MPa) give fourteen unknowns and fourteen equations. Oil and water always use the liquid-type EOS root; the vapor descendant uses the vapor-type root, including the unique-root continuation region. Starting pressures are guesses, not targets. No AAD is optimized. Temperature continuation preserves the nontrivial composition branch and rejects near-trivial collapse; rejected steps are recorded rather than certified as a critical point.

Independently, start the LL branch at 25 MPa for each temperature, decrease pressure, and locate the zero of its oil-phase constrained-composition Hessian minimum eigenvalue. This is a local composition spinodal, NOT a stable LL coexistence boundary. The central Hessian differentiates mu_i-mu_4 with x_j+=h and x_4-=h on the same liquid root. Curvature calculations are repeated at half h.

## Main result: intersection is NOT gas/oil critical merging

The retained vapor-onset locus crosses the independently computed oil spinodal near:

T = 589.045834 K, P = 11.615358 MPa.

Independent pressure-down spinodal refinement at the same temperature gives 11.6153582677 MPa; direct onset/curvature intersection gives 11.6153582747 MPa. Oil-parent composition L1 difference is 1.81e-10. Thus these are the same parent LL state within numerical error, not merely crossing interpolated plot lines.

At that intersection:

- oil water mole fraction: 0.719402827;
- vapor-descendant water mole fraction: 0.941380606;
- vapor/oil composition L1 distance: 0.443955558;
- vapor/oil MOLAR-density ratio: 0.705708420;
- oil minimum curvature approximately zero (2.46e-9; half-h -1.78e-8);
- vapor minimum curvature remains positive, approximately 6.687.

The finite composition/density separations rule out interpreting this intersection as gas/oil critical coalescence. It is a loss of local stability of the oil parent along a separately identified incipient-vapor stationary locus. It is not a certification of another liquid-liquid critical endpoint: higher-order and global phase-split tests are not performed here.

## Beyond the intersection the onset equations still have solutions, but the parent is unstable

| T K | retained vapor-onset P MPa | independent oil-spinodal P MPa | vapor/oil L1 | oil curvature on vapor locus |
|---:|---:|---:|---:|---:|
|583|10.655218|10.042292|0.507992|+0.737151|
|587|11.282085|11.078592|0.466190|+0.493595|
|589|11.607793|11.603287|0.444461|+0.018128|
|590|11.773905|11.867123|about 0.433|negative|
|593|12.286166|12.664090|0.398685|-4.180732|
|597|13.006239|13.737882|0.347247|-12.171849|
|600|13.585134|14.550107|0.291615|-18.993538|

The earlier high-temperature boolean-TPD boundaries sit close to the separate oil curvature locus: at 593 K, old 12.663320 versus spinodal 12.664090 MPa; at 603.6 K, old 15.530184 versus spinodal 15.530664 MPa; at 613.2 K, old 18.162182 versus spinodal 18.158925 MPa. These are not asserted to coincide exactly. They must not be spliced onto the 583 K vapor descendant and called one validated WLV-WL curve.

## A temperature fold, not a demonstrated critical merger

Ordinary prescribed-T continuation stalls around 600.2356 K / 13.6334 MPa. Adding temperature as an unknown and a pseudo-arclength equation traverses the turning region; 180 corrected arc points were recorded. The tangent temperature component changes sign and the curve returns to lower T.

Near the sampled temperature maximum, vapor/oil L1 remains about 0.2761 and the MOLAR-density ratio about 0.8923. The vapor-descendant curvature crosses zero, while the oil is already strongly locally unstable. This is an observed fold of stationary onset solutions, not gas/oil critical merging. No claim is made about the untraced remainder of the folded branch.

## Numerical checks and an important global-stability limitation

The 0.5 K and 0.25 K runs have forty common forward temperatures. Maximum pressure difference is 1.08e-11 MPa and maximum phase-composition difference is 3.05e-12. Maximum fourteen-equation residual over forward/reverse records is 1.12e-10. Reverse continuation reaches the 583 K anchor. Half-h oil curvature changes by at most 1.53e-7 along the onset runs and 1.26e-7 on the independent spinodal runs. Printed digits describe numerical consistency, not experimental uncertainty or physical-model accuracy.

Production unrestricted flash was also evaluated at onset P +/- 0.020 MPa for 583, 586, 588, 589, 589.1, 590 and 593 K. At 583/586/588/589 K the low side gives OGW and the high side OW, with positive minimum active-phase composition curvature and closed material balance in the tested states. This is supporting local/numerical evidence, not a proof of global minimization across arbitrary phase sets.

At 590 and 593 K on the low side the existing code can return three phases with missing_phase_stable=1 but NEGATIVE minimum active-phase composition curvature (approximately -0.07391 and -3.02189). Therefore missing-phase stability over unused O/G/W slots does not certify resistance to splitting an already active phase. These results must NOT be accepted as globally stable simply because all three slots are occupied.

Current interpretation:

- identified low-T vapor-onset segment: traced and locally checked;
- onset/spinodal intersection: numerically located, no gas/oil merger;
- post-intersection stationary onset and folded continuation: locally unstable, not a physical equilibrium boundary;
- actual stable high-T phase topology, including possible additional liquid splitting: unresolved;
- FIGURE7_REPRODUCTION = BLOCKED;
- FORMAL_FLOW = BLOCKED.

## Reproduce

From repository root:

```sh
flags='-std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -Iad/include -Icommon/include -Iindices/include -Imodels/include -Icase/include -Itest/include -Itools/include'
g++ $flags tools/example/cpa_branch_locus_audit/main.cpp -o /tmp/locus
g++ $flags tools/example/cpa_branch_locus_audit/refine.cpp -o /tmp/refine
/tmp/locus /tmp/locus_0p5 0.5
/tmp/locus /tmp/locus_0p25 0.25
/tmp/refine /tmp/locus_0p5/vapor_onset.csv /tmp/locus_refined
```

A successful diagnostic process means computation completed, not that every curve is stable or that a physical acceptance gate passed. Inspect event/status records and active curvatures. Generated CSVs and logs remain artifacts, not production parameter inputs.
