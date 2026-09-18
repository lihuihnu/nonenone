# H02 Run 30: completed trajectory, not accepted quantitative physics

Audit date: 2026-09-18.

## Provenance

- Repository: lihuihnu/nonenone, PR #2.
- Workflow run: 35355410117; job: 105633438205.
- Branch head used by the run: 9480a3b869b7409e41eeebed9c7d377a952c6224.
- Actual PR checkout recorded in artifact: de884442ef0ed57f220d1d009a01963ef28e58fa.
- Artifact: scw-h02-native-2d, ID 10551683103.
- Downloaded archive SHA-256: f1642429fe661214f85b495de87d83b628712f0eb15f0183b42519c6234a8ce8.
- Raw outputs and generated figures remain artifacts, not committed source.

The workflow finished successfully. A/B/C each completed 3030 accepted internal steps, zero rejected steps, and 6060 seconds of model time. All three final PVI values are approximately 2.02. The complete runs use a 20x8 verification grid; the 60x20 runs in this artifact remain two-step smoke tests. Injector control remained reservoir-total-rate without BHP switching.

## Corrected matched-PVI results

The old first-at-or-above selector used 2.02-PVI rows as the 2-PVI summary because the recorded 6000-second PVI was slightly below 2.0. The corrected analysis uses bracketing linear interpolation of actual PVI and records the endpoints and interpolation fraction. It neither extrapolates nor substitutes the last row when a target has not been reached. These are reprocessed native outputs, not new flow simulations.

| PVI | Mode | Heavy RF (%) | Heavy production (g) | Injector-producer BHP difference (kPa) | Producer Heavy-flux-weighted viscosity (mPa s) |
|---|---|---:|---:|---:|---:|
| 1 | A | 39.742312 | 22.897044 | 3.570723 | 2.000000 |
| 1 | B | 47.603533 | 27.426190 | 3.755759 | 2.000000 |
| 1 | C | 48.567481 | 27.981557 | 3.501243 | 1.732452 |
| 2 | A | 48.370782 | 27.868231 | 2.461032 | 2.000000 |
| 2 | B | 55.644744 | 32.059034 | 2.507163 | 2.000000 |
| 2 | C | 56.985567 | 32.831533 | 2.347616 | 1.732450 |

At 2 PVI, B-A is +7.273962 percentage points, C-B is +1.340824 points, and C-A is +8.614786 points. C-B reduces the reported producer viscosity by 13.377469% and the BHP difference by 6.363678%.

Viscosity here is reconstructed from the existing assumed H02 constitutive law and single-perforation phase mass fluxes. It is not a measured mixture viscosity and is not a domain-average viscosity. The C log-mixing law itself prescribes a decreasing viscosity as water mass fraction increases; this run measures the flow response to that assumption, not independent evidence that real Heavy must exhibit this reduction.

## Carrier audit

At 2 PVI:

| Mode | Oil-slot Heavy (g) | Gas-slot Heavy (g) | Water-slot Heavy (micrograms) | Water-slot share of total Heavy production |
|---|---:|---:|---:|---:|
| A | 27.868231 | 0 | 0 | 0 |
| B | 32.059033 | 0 | 0.719144 | 2.243186e-8 |
| C | 32.831533 | 0 | 0.711361 | 2.166700e-8 |

The phase-component cumulative sum matches the total ledger to below 2e-16 kg in the recomputation. A has zero Heavy through Water. B/C have zero gas saturation throughout the recorded reservoir diagnostics. At the 1- and approximately 2-PVI producer snapshots, the active mask is Oil+Water and the Water phase is essentially pure water. Thus the reported endpoints do not misclassify a mixed single phase as water-carried Heavy.

The Water-slot Heavy trace is thousands of times smaller than the maximum global Heavy mass defect. It is not a resolved quantitative extraction signal. The result does NOT establish substantial dissolved-Heavy carriage by SCW. The B-A response combines exchange, phase volumes/densities, saturations and flow formulation differences. It must not be called a pure extraction contribution. C-B is the response to the specified viscosity feedback protocol, still subject to numerical and physical validation.

## Strict mass audit: FAIL for B and C

The previous analyzer checked only final PVI, phase-ledger addition and zero Water-carried Heavy for A. It did not gate global mass balance. In addition, the common ledger's existing relativeError uses a scale with a 1-kg floor, inappropriate as the sole relative-error metric for this approximately 0.0576-kg initial Heavy inventory.

The revised postprocessor recomputes:

    e_i(t) = abs(M_i(t) - M_i(0) - Min_i(t) + Mout_i(t))
             / (M_i(0) + Min_i(t))

For Heavy this is initial-inventory-normalized error. Initially empty water is normalized by the actual cumulative injected water mass. Empty/negligible initial states are handled explicitly. All recorded output checkpoints are tested, not just the final row. This is not an every-internal-step inventory audit because those inventory snapshots were not saved.

| Mode | Maximum Heavy relative error | Maximum water relative error | Threshold | Verdict |
|---|---:|---:|---:|---|
| A | 4.151056e-9 | 2.409860e-8 | 1e-6 | PASS |
| B | 3.869223e-5 | 1.339496e-3 | 1e-6 | FAIL |
| C | 3.066861e-5 | 8.314634e-4 | 1e-6 | FAIL |

B/C maximum Heavy absolute defects are approximately 2.229205e-6 and 1.766934e-6 kg. The largest water relative errors occur at the first 60-second output checkpoint. The maximum Heavy relative errors occur at output steps 8 (B) and 25 (C). A small phase-ledger addition error does not certify global inventory conservation.

## Postprocessor changes and tests

- Reject insufficient PVI coverage rather than accepting 1.99 PVI as 2 PVI.
- Use bracketed actual-PVI interpolation rather than the first later row.
- Reject missing required phase columns and nonfinite required data instead of filling with zero.
- Preserve undefined/no-flow viscosity instead of forward-filling through gaps.
- Require one producer perforation for constitutive reconstruction and reject supported reverse-flow violations.
- Gate global inventory closure independently of phase-ledger closure.
- Keep all partial results labeled CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION.

Eight local regression tests passed, including re-auditing this exact artifact. With Run 30 inputs, the corrected analyzer writes the result tables and returns exit code 2 for the failed numerical mass gate. A future CI failure on the same trajectories is intentional evidence enforcement, not a claim that the native solver failed to run.

## Next bounded work

Do not tune EOS/BIPs or viscosity parameters. First diagnose the mass defect at the earliest saved 60-second checkpoint, separating nonlinear residual, post-solve phase-state update and inventory/source accounting. Run identical 20x8 cases with tightened numerical tolerances and a halved internal dt only as convergence controls. Then compare with complete 60x20 trajectories at matched actual PVI, with well mapping and Peaceman WI checked. No grid/time convergence is established by Run 30.

Current state:

    FULL_TRAJECTORY_COVERAGE = PASS
    PHASE_COMPONENT_LEDGER = PASS
    GLOBAL_MASS_NUMERICAL_GATE = FAIL
    GRID_TIME_CONVERGENCE = NOT_ESTABLISHED
    REAL_HEAVY_PHYSICAL_VALIDATION = NOT_VALIDATED

This audit changes only analysis and evidence handling. It does not modify production EOS, Heavy identity, kij=0.30, viscosity closures, geometry or solver tolerances, and it does not restart CPA/Figure-7 development.
