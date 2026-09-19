# Incipient identity audit: Gas slot does not prove vapor identity

## Scope and provenance

No Heavy/Athabasca parameter, EOS equation, production default, tolerance, or existing acceptance gate is changed by this diagnostic.

Audited physics revision: `07f9aaffe139e573bdf00637cf96b7c4ae8ec43e`.
Pointwise CI: run `35329819420` (run 35), artifact `10540751968`, ZIP SHA256 `e4c52b4546f61459f13f36ca3b3b02410e19f53064ee42e53d61cf1d80c2edca`.
Exact source: root-policy CI run `35329819474`, artifact `10539814168`, PR checkout `558133a79c74174d4e6773ec4fe6af1b0b7f83c1`.
High-T supplementary CI: run `35329819379`, artifact `10541006313`.
The added C++ diagnostic is commit `39e0af4fafce7b7dd8ca37acf683bb210b453947`. Its exact committed source was downloaded from artifact `10541515337` and compiled locally with C++17, O2, Wall/Wextra/Wpedantic/Werror. Source-file SHA256: `adfc5e93ad27001eeb21845508a97d8eb0688366681619fc72c89be2a9e55c8b`.

## Latest completed pointwise CI: path agreement 6/6, full gate 3/6

| T K | continued boundary MPa | independent/continued gap MPa | full pointwise gate |
|---:|---:|---:|---|
|522.9|4.158889119118|2.561160172831e-10|PASS|
|573.3|9.257909639296|1.629807400150e-10|PASS|
|583.0|10.65521773219|3.329478914793e-9|PASS|
|593.0|12.66332034383|2.165318591096e-9|FAIL: low-side global flash|
|603.6|15.53018432267|2.537854371099e-9|FAIL: low-side global flash|
|613.2|18.16218237060|5.238689482212e-9|FAIL: cold/continued OW branch|

All high-pressure-side global states are certified OW in this existing test. The existing full gate remains FAIL; no new pass criterion is substituted.

## New local experiment

At each of the six continued boundaries, evaluate P=Pb-offset for offsets -0.020, 0.0002, 0.001, 0.005 and 0.020 MPa. Continue the same fixed-liquid-root OW branch from 30 MPa. Recompute the candidate on selected/vapor/liquid roots, its composition distance from both active phases, molar-density ratios, TPD and stationary fugacity residual.

Independently form the oil-branch composition Hessian by differentiating mu_i-mu_4 under x_j+=h, x_4-=h on the SAME liquid root. Use h=1e-5*min(x_j,x_4), then repeat with half h. This is differentiation of production EOS chemical potentials, not a replacement thermodynamic model.

Both exact-source local runs completed 30/30 diagnostic states (60 evaluations across two finite-difference resolutions). Completion is NOT vapor-identity acceptance. Largest minimum-eigenvalue change after halving h: approximately 1.15e-7. Numerical Hessian asymmetry is at most approximately 1.31e-7.

At P=Pb-0.020 MPa:

|T K|candidate xH2O|oil xH2O|candidate/oil molar-density ratio|candidate-oil composition L1|oil Hessian minimum eigenvalue|
|---:|---:|---:|---:|---:|---:|
|522.9|0.969073|0.469026|0.309338|1.000094|+0.155611|
|573.3|0.952201|0.650903|0.585831|0.602596|+0.723536|
|583.0|0.946600|0.692564|0.653995|0.508072|+0.726342|
|593.0|0.733728|0.732468|0.998720|0.011169|-0.124023|
|603.6|0.766341|0.766456|0.996255|0.008578|-0.280874|
|613.2|0.795510|0.797109|0.992045|0.010451|-0.284343|

The density ratios above are MOLAR, not mass-density ratios. Eigenvalue magnitudes depend on the chosen independent composition coordinates; the negative sign is the relevant local-curvature evidence.

At 593 K the high-pressure-side oil Hessian minimum is +0.104488; it becomes -0.005525 at Pb-0.0002 MPa. At 603.6 K the corresponding values are +0.251404 and -0.009053. At 613.2 K it remains positive (+0.051287) at Pb-0.0002 and turns negative by Pb-0.005. Thus not all three high-T onsets should be called the exact same spinodal crossing.

## Interpretation and limitations

The first three states have an incipient phase distinct in both composition and molar density from the parent oil. The high-T candidates instead approach the parent oil in composition and density, and the restricted oil branch exhibits negative local composition curvature below the reported boundary. This is evidence of oil-like compositional splitting/critical behavior contaminating an unconditional WLV-WL interpretation; it is NOT proof that every high-T state is LLL, nor proof that vapor is absent.

A unique density root alone neither proves nor disproves vapor identity. In particular, candidate_Z==vapor_Z==liquid_Z cannot identify a vapor branch. A Gas slot, tiny negative TPD, and a converged three-slot solution do not replace continuation of a physically identified branch. No arbitrary density-ratio cutoff is introduced as a physical acceptance tolerance.

Next discriminating test: continue the well-separated 583 K incipient vapor in temperature and pressure while retaining root/composition identity, separately track the oil-like composition-instability branch, and check active-phase curvature/global Gibbs stability before deciding which locus is the actual WLV-WL boundary. Do not compensate for this ambiguity by fitting to aggregate AAD or merely increasing recovery iterations.

Status: `UNRESOLVED_HIGH_T_NEAR_OIL_COMPOSITION_SPLIT`; full Figure-7 reproduction and flow promotion remain BLOCKED.

## Reproduce

From repository root, using the downloaded run-35 figure7_pointwise.csv as INPUT:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror \
  -Iad/include -Icommon/include -Iindices/include -Imodels/include \
  -Icase/include -Itest/include -Itools/include \
  tools/example/cpa_incipient_identity_audit/main.cpp -o /tmp/identity_audit
/tmp/identity_audit "$INPUT" /tmp/identity_h 1e-5
/tmp/identity_audit "$INPUT" /tmp/identity_h2 5e-6
python tools/example/cpa_incipient_identity_audit/analyze.py \
  /tmp/identity_h /tmp/identity_h2 /tmp/identity_analysis
```

Generated CSVs, logs, binaries and figures are artifacts, not committed source. The new diagnostic was executed locally; the provenance above distinguishes that execution from the pre-existing GitHub pointwise CI.
