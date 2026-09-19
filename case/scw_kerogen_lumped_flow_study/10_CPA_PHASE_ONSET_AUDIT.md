# 10 — CPA phase-onset continuation audit

## Scope

This audit isolates the two CPA zero-dimensional failures originally observed at

- `T = 653.15 K (380 °C)`;
- BASE oil composition;
- `z_H2O = 0.7610625`;
- `p = 26.00` and `26.25 MPa`.

All CPA parameters are frozen. No `a0`, `b`, `c1`, association parameter, cross-association parameter, or `kij` is changed.

The diagnostic executable is:

`tools/example/scw_cpa_onset_audit/main.cpp`

It uses the production CPA EOS, production flash, public stability test, restricted active-set solver, canonical phase roles and thermodynamic profiler.

## Hypotheses tested

Five explanations were tested independently:

1. pressure-continuation resolution / seed reuse near phase onset;
2. active-set switching or pair initialization;
3. Oil/Water role canonicalization;
4. CPA association-site or density-root numerical failure;
5. absence of a physically valid solution under the current CPA parameter/model form.

The final classification is:

| hypothesis | result |
|---|---|
| CPA association / density-root numerical failure | **not supported** |
| Oil/Water role canonicalization failure | **not supported** |
| unseeded active-set / O+W pair initialization | **supported** |
| continuation / seed reuse deficiency | **supported** |
| CPA parameter/model-form pathology required to explain these failures | **not required** |

This classification concerns these two numerical failures only. It does not promote the current CPA parameters to experimentally calibrated status.

## Evidence

### 1. The Water-only state is not thermodynamically stable

At both failed states, a Water-only restricted result exists numerically, but a fresh stability test rejects it.

At 26.00 MPa:

- Oil trial sum ≈ `1.000142`;
- Gas trial sum ≈ `1.000142`;
- both missing Oil and missing Gas are flagged unstable.

At 26.25 MPa:

- Oil trial sum ≈ `1.000325559453018`;
- Gas trial sum ≈ `1.000325559453019`;
- both missing Oil and missing Gas are flagged unstable.

Therefore the failure is not evidence for a stable single aqueous phase.

### 2. The Oil and Gas instability directions are numerically the same direction

The missing Oil/Gas trial strengths are equal to roundoff, and their incipient compositions are also coincident.

The incipient Oil/Gas composition L1 distances are approximately:

- 26.00 MPa: `1.92e-14`;
- 26.25 MPa: `2.58e-14`.

Thus the Water-only stability calculation is reporting the same nonaqueous bifurcation through two public nonaqueous roles in a CPA one-density-root region. Treating the two labels as independent phases causes an artificial one-phase -> three-phase active-set jump.

### 3. CPA density roots and association are healthy

Direct CPA root probes at 25.75, 26.00, 26.25 and 26.50 MPa successfully evaluate Oil, Gas and Water roles.

At the two failed states:

- all requested roots are finite and positive;
- the 4C-water association calculation stays on the analytic/reduced fast path;
- iterative association calls are zero;
- association iterations are zero.

This excludes association-site iteration or inability to bracket a CPA density root as the cause of these failures.

### 4. Canonical Oil/Water roles are healthy

Using a nearby physical O+W state as a seed, then deliberately swapping the Oil and Water seed compositions, both restricted solves return the same canonical tie-line.

The final composition L1 discrepancy is of order `1e-11` or smaller. Therefore role canonicalization is not the failure mechanism.

### 5. A stable O+W equilibrium exists with the frozen parameters

Using a nearby converged O+W composition only as an initial iterate, the same frozen CPA model converges at both target states.

After the production fix, unrestricted flash now recovers these states directly.

At 26.00 MPa:

- phase code = `5` (O+W);
- `beta_o ≈ 0.986800`;
- `beta_w ≈ 0.013200`;
- `x_H2O^o ≈ 0.759063`;
- `x_H2O^w ≈ 0.910556`;
- component material closure ≈ `1.1e-16`;
- max active-phase log-fugacity spread ≈ `1.8e-14`;
- final stability certificate = PASS.

At 26.25 MPa:

- phase code = `5` (O+W);
- `beta_o ≈ 0.972019`;
- `beta_w ≈ 0.027981`;
- `x_H2O^o ≈ 0.756695`;
- `x_H2O^w ≈ 0.912780`;
- component material closure ≈ `1.1e-16`;
- max active-phase log-fugacity spread ≈ `2.0e-13`;
- final stability certificate = PASS.

A physically certified solution therefore existed before any parameter change. The original non-convergence was a path/active-set problem.

## Production fix

The production active-set logic now handles a narrowly defined CPA case:

1. active set is Water-only;
2. missing Oil and Gas are both unstable;
3. their trial strengths agree within the stability tolerance scale;
4. their incipient compositions are coincident within the phase-composition tolerance.

For that case the solver first stages the liquid-root O+W branch using the Oil incipient composition.

- If an O+W two-phase basin exists, it is reused and stability is recomputed. A genuinely distinct Gas instability can then be released on the next outer iteration.
- If the O+W branch collapses back to one phase, the solver restores the original simultaneous Oil+Gas release so a genuine O+G+W onset is not suppressed.

This is a numerical active-set change only. The thermodynamic model and all CPA parameters remain unchanged.

## Regression gates

The fix passes all of the following:

- exact 26.00/26.25 MPa failure-state unit regression;
- existing deterministic PR/SW/CPA thermodynamic failure-mode sweep;
- public CO2–DME–H2O CPA VLLE benchmark:
  - all four published three-phase states recovered;
  - composition MAE remains approximately `0.01387`;
  - maximum absolute composition error remains approximately `0.04156`;
- full SCW kerogen 0D acceptance:
  - PR registered scan: PASS;
  - CPA registered scan: PASS;
  - PR P-T map health: PASS;
  - CPA P-T map health: PASS;
  - PR dense composition paths: PASS;
  - CPA dense composition paths: PASS;
  - cross-EOS initial-state gate: PASS.

The strict runtime result is now:

`ZERO_D_PVT_GATE_PASS`.

## Interpretation

The two original CPA failure states are classified as an **active-set / initialization and continuation robustness defect**, not as evidence that the frozen CPA parameter set lacks an equilibrium solution.

This result does **not** establish that the CPA parameter set is experimentally correct. The independent CPA calibration, density validation and viscosity validation gates remain authoritative.

Therefore passing this numerical 0D gate removes one pre-flow numerical blocker, but reservoir flow remains blocked until the remaining experimental/property gates pass.
