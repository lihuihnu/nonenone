# SCW binary zero-flow calibration

This offline executable calibrates the water–squalane binary without entering
the reservoir flow loop. It reuses the production PR EOS, three-phase flash,
volume-translation path and LBC viscosity implementation.

The calculation is deliberately split into three ordered regressions:

1. fit the squalane component volume translation to pure-squalane density;
2. hold that translation fixed and fit the effective LBC critical volume to
   pure-squalane viscosity;
3. fit `k_H2O-squalane(T) = k_ref + b(1/T - 1/653.2 K)` directly to the
   production-flash coexistence compositions.

Alternating temperatures are used for training and validation in the
pure-squalane table. The two LLE calibration rows are 637.2 K/26.88 MPa and
653.2 K/27.74 MPa; the other two rows are held out. Reported experimental
accuracy bands are deterministic acceptance limits, not standard deviations.

The pure-water check independently evaluates the industrial form of
IAPWS-2008 (`mu2=1`) using production IF97 density, and compares it with the
locked McBride–Wright closure. The IAPWS evaluator is verified against the
official 298.15 K, 998 kg/m3 sample point before any output is created.

Run from the repository root in WSL/Linux:

```bash
make -C tools run-scw-binary-calibration CXX=clang++
```

The target refuses to overwrite an existing `results/` directory. Generated
results are ignored by Git. See `docs/SCW_BINARY_ZERO_FLOW_CALIBRATION.md` for
the accepted/rejected parameter decisions and applicability limits.

