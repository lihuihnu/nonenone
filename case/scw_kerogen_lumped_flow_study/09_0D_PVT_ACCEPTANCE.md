# 09 — Complete zero-dimensional PVT acceptance before flow

## Hard rule

No reservoir flow comparison is allowed to decide which EOS is physically reasonable.

PR and CPA must first be exercised **independently** on the same zero-dimensional PVT state space. Only after both models return physically self-consistent initial states, and the pre-existing phase/density/viscosity gates are satisfied, may a flow comparison be called fair.

The executable is:

`tools/example/scw_kerogen_0d_pvt_acceptance/main.cpp`

It reuses the production EOS, production three-phase flash, public stability test, and phase-diagram utilities. No second flash or second EOS is implemented.

## Registered target temperatures and pressures

The exact acceptance temperatures are:

- 360 °C = 633.15 K;
- 374 °C = 647.15 K;
- 380 °C = 653.15 K.

The target pressure grid is:

- 25, 26, 27, 28, 29, 30 MPa.

A denser envelope diagnostic uses 628.15–658.15 K and 20–35 MPa so that all three target temperatures and the 25–30 MPa region lie inside the map.

## Overall-composition scan

The central oil composition is the four-lump recovered-oil mole distribution from the current characterization:

- Gasoline = 0.036330;
- Diesel = 0.406770;
- Middle = 0.325421;
- Heavy = 0.231479.

For each oil-composition family, the H2O overall mole fraction is scanned over:

`0.01, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 0.97, 0.995`.

Three deterministic oil families are registered:

1. `BASE` — measured characterized oil mole ratio;
2. `LIGHT_ENRICHED` — deterministic light-side sensitivity;
3. `HEAVY_ENRICHED` — deterministic heavy-side sensitivity.

The latter two are **not uncertainty bounds and not experimental compositions**. They are pre-registered stress tests intended to reveal phase-role switching, heavy-root problems and composition-local flash pathologies.

The machine-readable grid is `pvt_acceptance/composition_scan.csv`.

## Initial-state anchor

The fairness anchor is the `BASE` oil distribution at

`z_H2O = 0.20`

and

`p = 25 MPa`.

It is evaluated independently at all three target temperatures. The corresponding overall oil composition is kept identical between PR and CPA. Each EOS is allowed to predict its own equilibrium phase count and phase fractions.

**The two EOS are not required to predict the same phase count.** A phase-count difference is a scientific model result, not an automatic failure. What is required is that each predicted equilibrium is internally physical.

## State-level physical consistency

For every registered T–P–z state the acceptance harness requires:

1. production flash converges;
2. at least one phase is active;
3. the final active set passes a fresh production `stabilityTest()`;
4. component material reconstruction error is <= `1e-8`;
5. maximum active-phase log-fugacity spread is <= `1e-6`;
6. every active phase composition is finite, non-negative and normalized;
7. every active phase has positive finite Z, molar density, mass density and LBC viscosity;
8. if a Water-role phase exists, it is the most water-rich active phase.

The role convention is the production canonical `Oil / Gas / Water` convention.

## Stability output

`state_scan.csv` stores:

- final phase code and phase count;
- stability-valid / stability-stable flags;
- missing-phase unstable flags;
- TPD-style trial sums for Oil/Gas/Water candidates;
- mass closure;
- fugacity closure;
- role and property checks;
- final per-state pass/fail.

This means a flash result is not accepted merely because Newton converged.

## Phase composition, density, viscosity and phase role

`phase_properties.csv` writes one row for each canonical role at every state:

- active/inactive flag;
- phase mole fraction;
- saturation;
- compressibility;
- molar density;
- mass density;
- LBC viscosity;
- IAPWS-2008 water viscosity when the Water-role composition remains inside the explicit <=2 mol% solute domain;
- all five phase mole fractions.

The IAPWS value is an auxiliary water-rich reference. It does not replace the common transport closure or silently extrapolate into hydrocarbon-rich water-role states.

## Phase-envelope products

For PR and CPA separately, and for BASE/LIGHT_ENRICHED/HEAVY_ENRICHED at `z_H2O=0.20`, the harness writes:

- unrestricted O/G/W P–T maps;
- refined oil/gas/water phase-onset boundaries;
- a clearly labelled restricted O/G bubble/dew **projection**.

The O/G projection is not substituted for the full aqueous phase envelope.

The harness also writes dense pressure-composition maps from oil-rich to water-rich composition at 360, 374 and 380 °C over 25–30 MPa.

## Acceptance levels

There are two distinct runtime gates.

### Model-domain scan gate

Every registered target T–P–z state for a backend must pass the state-level physical-consistency checks.

This is intentionally stronger than testing only the initial point because a compositional flow trajectory can enter enriched/depleted states.

### Cross-EOS initial-state gate

At 360/374/380 °C and 25 MPa, the `BASE, z_H2O=0.20` state must pass independently for both PR and CPA.

The output `cross_eos_initial_state.csv` reports phase-count agreement only as a diagnostic.

The runtime zero-dimensional gate is PASS only if:

`PR scan PASS && CPA scan PASS && both-EOS initial-state PASS`.

## Reservoir-entry dependency

A zero-dimensional structural PVT pass is necessary but not sufficient for reservoir flow.

The final flow-comparison readiness also depends on:

- accepted H2O–lump PR binary calibration;
- accepted CPA binary/association calibration;
- density validation;
- viscosity validation.

The committed dependency table is `pvt_acceptance/flow_entry_gate.csv`.

Thus a numerically clean 0D scan cannot override missing experimental calibration, and a good flow result cannot override a failed 0D state.

## Current runtime result

The GitHub-run production-kernel acceptance scan has now been executed.

| gate | PR | CPA |
|---|---:|---:|
| registered stability/PVT states | 594 / 594 PASS | 594 / 594 PASS |
| target-window P-T envelope-map flash states | 5673 / 5673 converged | 5673 / 5673 converged |
| dense pressure-composition states | 5103 / 5103 converged | 5101 / 5103 converged |
| BASE initial anchors at 25 MPa | 3 / 3 PASS | 3 / 3 PASS |

All three BASE initial anchors (360/374/380 °C, `z_H2O=0.20`, 25 MPa) are internally self-consistent in both EOS. Both currently return a single Oil-role phase.

The complete 0D gate is nevertheless **BLOCKED** because CPA has two local non-converged states on the dense 380 °C composition path:

- path fraction `0.7625`, corresponding to `z_H2O ≈ 0.7610625`;
- 26.00 MPa;
- 26.25 MPa.

Neighbouring CPA states show that these failures lie at the Water-only -> Oil+Water onset region: the same composition is Water-only below the narrow interval and O+W above it. This is exactly the kind of phase-boundary robustness issue that a compositional flow calculation may encounter, so it is not waived and is not repaired by BIP tuning.

The PR dense path has no corresponding failure.

The full BASE `z_H2O=0.20` P-T maps over 628.15–658.15 K and 20–35 MPa contain no phase-onset boundary for either EOS; both remain single Oil-role over that local initial-composition window. This does not imply the full composition space is single phase.

The coarse registered composition scan also shows a substantive EOS topology difference: PR returns one- or two-phase O/W states in the registered domain, whereas CPA additionally produces three-phase O/G/W states in high-water, high-temperature / lower-pressure parts of the scan. This difference is retained as a model result rather than ranked as better or worse.

At the initial anchors, CPA currently predicts mass density about 9.8–9.9% below PR and LBC viscosity about 25% below PR. Because the independent density and viscosity validation gates remain blocked, these differences are **not** yet admissible as flow-model evidence.

Therefore:

- `CROSS_EOS_INITIAL_STATE = PASS`;
- `PR_DENSE_COMPOSITION_PATH = PASS`;
- `CPA_DENSE_COMPOSITION_PATH = FAIL`;
- `ZERO_D_PVT_ACCEPTANCE = BLOCKED`;
- reservoir flow remains frozen.

The next thermodynamic task is to diagnose the two CPA phase-onset failures with production-flash continuation / boundary refinement while keeping the current CPA physical and BIP parameters fixed. Only after numerical robustness is separated from parameter/model-form error should any parameter calibration resume.
