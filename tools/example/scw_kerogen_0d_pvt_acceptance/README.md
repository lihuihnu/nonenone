# SCW kerogen 0D PVT acceptance

This zero-flow executable is the final thermodynamic preflight before any PR-vs-CPA reservoir comparison.

It reads the committed case parameter files rather than maintaining a second hidden parameter set:

- PR pure parameters and PR BIP matrix;
- CPA pure/association parameters and CPA BIP matrix;
- the registered composition scan.

It then runs the **production** `CubicEquationOfState`, `CubicThreePhaseFlash`, `stabilityTest` and `PhaseDiagramSampler`.

## Run

From repository root:

```bash
make -C tools run-scw-kerogen-0d-pvt-acceptance CXX=clang++
```

The normal run target executes the study and preserves a blocked scientific gate as a successful tool execution. To use it as a hard prerequisite:

```bash
make -C tools require-scw-kerogen-0d-pvt-acceptance CXX=clang++
```

The `require` target propagates exit code 3 when either EOS has a failed registered state or the cross-EOS initial-state gate is blocked.

Generated `results/` are ignored by Git.

## Outputs

Each backend gets:

- `state_scan.csv`
- `phase_properties.csv`
- `acceptance_summary.csv`
- full P-T maps and phase-onset envelopes for the three oil-composition families
- restricted O/G envelope projections
- pressure-composition maps at 360, 374 and 380 °C

The root result directory gets:

- `cross_eos_initial_state.csv`
- `zero_d_pvt_gate.csv`
- `zero_d_pvt_gate.txt`

See `case/scw_kerogen_lumped_flow_study/09_0D_PVT_ACCEPTANCE.md` for the scientific acceptance contract.
