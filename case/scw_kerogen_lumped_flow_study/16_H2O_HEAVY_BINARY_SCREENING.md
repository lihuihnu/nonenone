# 16 — H2O–Heavy binary screening before flow

## Scope

This stage is the deliberately minimal constant-temperature system proposed before adding Light/Middle components:

`H2O + OIL_HEAVY`

No kerogen cracking, upgrading reaction, coke formation, component generation or reaction kinetics is enabled. The purpose is to isolate the thermodynamic/transport consequences of the current single Heavy pseudo-component and to expose which missing measurements prevent a physical flow interpretation.

The screening is executed inside the existing five-component **production** kernel. Gasoline, Diesel and Middle have exactly zero total inventory. Heavy is `1-z_H2O`.

## Frozen Heavy identity used for this screening

The current `OIL_HEAVY` is not squalane. It is the provisional `>500 °C` tail pseudo-component from the experiment-driven characterization:

| property | current value |
|---|---:|
| MW | 660.132 g/mol |
| average Tb | 563.696 °C |
| SG | 0.942605 |
| Tc | 982.864 K |
| Pc | 0.808714 MPa |
| acentric factor | 1.215676 |
| Vc | 1793.890 cm3/mol |
| characterization status | `PROVISIONAL_HEAVY_TAIL_LOWER_BOUND_LIKE` |

The recovered ACS 2023 main article directly defines the simulated-distillation Heavy cut as `>500 °C`; the paired SI supplies the corresponding fractions. This is useful for pseudo-component identity but is not H2O–Heavy phase-equilibrium calibration.

## Interaction-model status

The current binary parameters remain frozen during this screening:

- PR H2O–Heavy: the current numerical prior is inherited from the **water+squalane benchmark** and is explicitly `HEAVY_BENCHMARK_ONLY`; it is not a fit to the generated Heavy pseudo-cut.
- CPA H2O–Heavy: `kij=0` is diagnostic only.
- Heavy self-association and H2O–Heavy cross-association are not asserted in the baseline because the required Heavy PNA/polar/resin/asphaltene evidence is not available.

No BIP or CPA association parameter was changed to obtain the numerical results below.

## Production-kernel screening design

Registered temperatures:

- 360 °C = 633.15 K;
- 374 °C = 647.15 K;
- 380 °C = 653.15 K.

Registered pressures:

`25, 26, 27, 28, 29, 30 MPa`.

Registered H2O mole fractions:

`0.01, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 0.97, 0.995`.

The production acceptance tool also runs:

- a full O/G/W P–T map at `z_H2O=0.20` over 628.15–658.15 K and 20–35 MPa;
- dense H2O–Heavy pressure–composition paths with 81 composition points × 21 pressure points × 3 target temperatures = **5103 states per EOS**.

## Lower-dimensional fugacity-closure correction

The first registered binary run correctly exposed an acceptance-harness defect rather than a parameter defect.

The flash itself converged, passed stability, material closure, role checks and finite-property checks, but the acceptance residual still required equal fugacity for Gasoline, Diesel and Middle even though their total inventories were exactly zero. In a lower-dimensional mixture, a zero-inventory component is governed by a complementarity condition; equality of its trace numerical fugacities across active phases is not an active equilibrium equation.

The acceptance check was therefore corrected to apply the `1e-6` log-fugacity equality tolerance only to components with `z_i > 1e-14`. Material closure and non-negative normalized phase compositions still include all five components. No EOS or BIP parameter was modified by this correction.

## Numerical result

The production-kernel structural screen now passes:

| numerical check | PR | CPA |
|---|---:|---:|
| registered H2O–Heavy states | 198 / 198 PASS | 198 / 198 PASS |
| dense H2O–Heavy P–z states | 5103 / 5103 converged | 5103 / 5103 converged |
| full H2O–Heavy target-window P–T map | no flash failures | no flash failures |
| maximum present-component log-fugacity spread in the registered binary states | ~1.25e-12 | ~1.10e-12 |

This establishes **numerical admissibility only**.

## Current phase-split prediction at 28 MPa

The dense composition grid brackets the first one-phase → two-phase transition as follows:

| T | PR onset bracket in overall zH2O | CPA onset bracket in overall zH2O |
|---|---:|---:|
| 360 °C | 0.613313–0.625625 | 0.440937–0.453250 |
| 374 °C | 0.650250–0.662563 | 0.502500–0.514813 |
| 380 °C | 0.662563–0.674875 | 0.539438–0.551750 |

These are grid brackets, not fitted experimental phase boundaries.

At 380 °C and 28 MPa the current two-phase tie-line-like screening result is:

| quantity | PR | CPA |
|---|---:|---:|
| H2O mole fraction in oil-rich phase | 0.668194 | 0.541401 |
| H2O mass fraction in oil-rich phase | 5.21 wt% | 3.12 wt% |
| Heavy mole fraction in water-rich phase | 3.73e-9 | 5.54e-10 |
| oil-rich density | 655.75 kg/m3 | 634.07 kg/m3 |
| oil-rich LBC viscosity | 0.1957 mPa s | 0.1825 mPa s |
| water-rich density | 364.46 kg/m3 | 570.86 kg/m3 |
| water-rich LBC viscosity | 0.0469 mPa s | 0.0848 mPa s |
| pure-water IAPWS viscosity reference at this T/P | 0.05935 mPa s | 0.05935 mPa s |

The oil-rich water mole fractions look large because the Heavy pseudo-component has MW ≈660 g/mol; on a mass basis the predicted dissolved-water levels are only about 3–5 wt%.

The PR/CPA onset difference is itself a warning: with unvalidated H2O–Heavy parameters, model choice materially changes the predicted phase split.

More importantly, both baselines put essentially zero `>500 °C` Heavy in the water-rich phase. This means the current single-Heavy screening model does **not** yet provide a defensible model of the experimentally reported SCW extraction/carrying mechanism. It must not be used to claim that SCW carries or fails to carry the actual generated Heavy fraction.

## Literature evidence hierarchy

### Direct equilibrium context — not target-window calibration

Matsui et al. (2014), DOI `10.1627/jpi.57.118`, measured water + Canadian bitumen equilibrium at 603–653 K and up to about 16–17 MPa. At 653 K the visual observations report formation of a new water-rich liquid near 17 MPa and disappearance of an observable vapor by about 20 MPa; at 25 MPa the view cell became too dark for a quantitative phase observation.

This is valuable qualitative high-temperature phase-topology context, but the quantitative equilibrium data do not reach the 25–30 MPa target window and the oil is not the generated `OIL_HEAVY` pseudo-cut.

Sato et al. (2018), DOI `10.1627/jpi.61.256`, measured water + atmospheric residue VLE at 603–643 K and 2.0–10.2 MPa. Their PR analysis showed substantially better liquid-water prediction when the residue characterization used detailed molecular-structure information rather than only average MW/SG. This supports treating Heavy characterization uncertainty as a first-class model issue.

### Flow/mechanism context — never a BIP fitting target

Zhao et al. (2018), DOI `10.1021/acs.energyfuels.7b03839`, reported SCW flooding of extra-heavy oil at 400 °C and 25 MPa and attributed the recovery mechanism to extraction of heavy-oil components into a water-rich phase and miscible flooding.

Zhao et al. (2020), DOI `10.1021/acs.energyfuels.9b03946`, further reported SCW injection as both heat carrier and organic solvent and discussed heavy-oil miscibility/upgrading.

These studies support the **existence of a plausible SCW extraction/transport mechanism**, but they do not provide the target Heavy pseudo-component equilibrium tie-lines required to regress this project's PR/CPA H2O–Heavy parameters.

The machine-readable source roles are in `binary_pr_calibration/h2o_heavy_evidence_manifest.csv`.

## Physical gate

The binary numerical screen does not promote flow. The following remain mandatory:

1. target-window H2O–actual-Heavy phase count / phase boundary / tie-line composition data near 360–380 °C and 25–30 MPa;
2. independent Heavy density data at the same window;
3. independent Heavy or reconstructed-oil viscosity data at the same window;
4. chemical evidence (PNA/aromatic/polar/resin/asphaltene/water-solubility) before adding Heavy self-association or H2O–Heavy cross-association in CPA;
5. hold-out rows not used for parameter regression.

Until these are available:

`H2O_HEAVY_PHYSICAL_VALIDATION = BLOCKED`

and

`FLOW_PROMOTION = BLOCKED`.

## Next numerical experiment after data acquisition

Once target-window Heavy data exist, fit PR and CPA independently to the calibration split, freeze parameters, and compare the same hold-out rows. Only then use the binary model in the laboratory slab to study:

- equilibrium oil-/water-rich viscosities;
- Heavy transfer into the water-rich phase;
- production history under SCW injection.

After that binary mechanism is credible, add one lighter oil pseudo-component to form `H2O + Light + Heavy` and test selective extraction/compositional shift at the producer.
