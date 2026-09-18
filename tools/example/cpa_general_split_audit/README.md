# Active-phase stability and variable-phase Gibbs splitting

## Scope

This independent diagnostic/research solver vetoes acceptance of unstable occupied phases and can represent multiple liquid-root records plus a vapor-root record. It does NOT replace the three-slot production flow state type. No OIL_HEAVY or Jia/Athabasca parameter, production EOS equation, production tolerance, or production default is changed. FIGURE7_REPRODUCTION and FORMAL_FLOW remain BLOCKED.

The EOS is evaluated by a C++ shared-library bridge that includes the exact frozen `cpa_athabasca_bitumen_water/main.cpp` factory. Python does not implement or approximate CPA. The bridge alone still reads the old O/G/W flash as an INITIAL GUESS. Subsequent minimization uses a variable-length inventory matrix, not those three slots.

Local source provenance: artifact 10541515337, source checkout 303bda33d2ef5e67e6a10c7300f25ae6250e626b (physics source 39e0af4fafce7b7dd8ca37acf683bb210b453947). The comparison to 78b3c2049f799f7fe93711fc60efaade4c5e5875 contains only subsequent diagnostics/documentation. Frozen factory Git blob: f23ad37e7673e756274e36762ce7b3f43aa1c3ad. Local compiler: GCC 14.2, C++17, O2, Wall/Wextra/Wpedantic/Werror; numerical environment: SciPy 1.17.0. Remote CI independently records its checkout and installed dependencies. Local numbers below must not be relabelled as remote CI results.

## Algorithm and acceptance boundary

1. For EVERY occupied phase, differentiate production chemical-potential differences at a fixed root on the normalized composition simplex. Repeat at half the finite-difference step. Negative minimum curvature exceeding the numerical uncertainty is a rejection, not a new absent-slot label. Near-zero curvature is separately recorded as unresolved.
2. For each active-phase chemical-potential reference, search liquid AND vapor TPD candidates, including root families already occupied. No missing-slot filter is used. A negative evaluated TPD is a witness even if an optimizer stops early; incomplete searches are counted and prevent global-stability acceptance.
3. A negative-curvature phase seeds two same-root daughters along its unstable eigendirection, with exactly conserved component inventory. Finite-distance negative TPD can also seed a new record by removing its inventory from the reference parent. Only a directly verified decrease in total Gibbs energy permits seed insertion.
4. For q phase records and N components, write n[k,i]=z[i]*softmax(u[:,i])[k]. This gives nonnegative inventory and exact component conservation at every iterate. Minimize sum(n[k,i]*mu[k,i]) with analytic allocation gradients. q is not fixed at three; the diagnostic resource ceiling is six records and never authorizes physical acceptance.
5. Mass-weighted optimization gradients can hide trace-component fugacity mismatch. Therefore a separate chemical-potential corrector enforces all (q-1)*N interphase equalities while retaining exact mass allocation. Reject an upward Gibbs change. Recompute curvature, material closure, fugacity equality and candidate TPD after relaxation.

The root family is attached to each record, not used as a unique phase name. Repeated liquid-root records are allowed. Root labels do not establish experimental gas/liquid identity in unique-root regions. A positive local Hessian plus a finite TPD search is not a proof of the global minimum.

## Actual completed local regression

Three CPA states use 44.1 wt% water / 55.9 wt% Athabasca bitumen. Pressures are the preceding independently continued vapor-onset pressure minus 0.020 MPa; they are not refitted here.

| T K | P MPa | old -> new record count | original minimum active curvature | final minimum active curvature | delta G/(nRT) | final max log-fugacity spread |
|---:|---:|---:|---:|---:|---:|---:|
|588|11.423873867918125|3 -> 3|+0.47654214|+0.47654214|0|1.14881e-10|
|590|11.753904875999448|3 -> 4|-0.07390917|+0.10505901|-3.05788547e-8|1.42109e-14|
|593|12.266165902108284|3 -> 4|-3.02189320|+1.01063841|-1.39094392e-5|4.97380e-14|

All three original states report legacy missing-phase stability. The new occupied-phase audit rejects the 590/593 K results despite that flag. Both yield four-record stationary candidates with three liquid-root records (two hydrocarbon-rich plus water-rich) and one vapor-root record. The lower total Gibbs energy at identical T/P/z is a direct witness that the original three-record state is not the global minimum.

Maximum final component material-closure error among these cases: 1.11023e-16. Reversing the unstable eigenvector reproduces the same daughter set up to permutation within 1e-6 composition, with total Gibbs differences 1.11023e-16 at 590 K and 4.44089e-16 at 593 K.

At 593 K the two hydrocarbon-rich daughter compositions differ especially in the Jia asphaltene pseudo-component: mole fractions approximately 0.01235483 and 0.04457869. These are NOT OIL_HEAVY characterization inputs. Their amounts are approximately 0.08483768 and 0.06944178 per normalized total feed mole; vapor-root and water-rich amounts are 0.00107623 and 0.84464431.

Analytic regular-solution regression also passes: ideal-mixture positive curvature; symmetric demixing negative curvature and lower-G mass-conserving split; and a locally convex but finite-distance TPD-unstable occupied liquid. This checks that local positive curvature is not incorrectly treated as sufficient stability.

## What remains unresolved

The initial BFGS stage stops with precision-loss status in the two physical split cases; this is explicitly retained. The independent chemical-potential corrector converges in 7 and 5 evaluations, and the final state is checked rather than trusting the BFGS success flag.

No materially negative final TPD witness was found in the sampled searches (minimum about -3.07e-15 at 590 K and -1.59e-15 at 593 K), but 10 and 9 optimizer searches respectively remain incomplete. Consequently both are explicitly UNRESOLVED_TPD_SEARCH, not globally certified equilibrium. The control also retains incomplete-search status. Do not relabel them PASS by discarding failed starts.

This proves that better, four-record, locally composition-stable stationary candidates exist; it does not prove the real experimental system has four phases. Complete global candidate coverage, mechanical/root admissibility, removal/merging of vanishing or duplicate records, and general phase-variable storage/transport integration remain separate tasks. No four-record state is silently truncated to the production three-slot layout.

## Reproduce

From repository root:

```sh
mkdir -p /tmp/active-split
python -m pip install numpy scipy
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -fPIC -shared \
  -Iad/include -Icommon/include -Iindices/include -Imodels/include \
  -Icase/include -Itest/include -Itools/include \
  tools/example/cpa_general_split_audit/eos_bridge.cpp -o /tmp/libsplit.so
python tools/example/cpa_general_split_audit/verify.py /tmp/libsplit.so /tmp/active-split
```

`ALGORITHM_REGRESSION_PASS` means the registered rejection/splitting/conservation tests pass. It does NOT change FIGURE7_REPRODUCTION=BLOCKED or FORMAL_FLOW=BLOCKED. Outputs include every phase composition, amount, fixed-root Z, chemical potential, curvature at both steps, optimizer status and incomplete TPD-search counts.
