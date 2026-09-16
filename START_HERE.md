# Start here

This is the unified, slimmed MPMC_SCW development repository. The authoritative
state is the Git `main` branch plus the executable tests; generated figures,
run results, external grids and machine-local dependencies are intentionally
not bundled.

## First checkout on any machine

1. Read `README.md`, `AGENTS.md` and `config/README.md`.
2. Copy the closest profile to `config/hpc.local.mk` and edit dependency paths.
3. Run `make print-config` and the appropriate `make doctor-*` command.
4. Run `make audit` and `make unit`.
5. Compile one representative case before starting long production runs.

For the H2O-CO2-nC10 benchmark:

```bash
make case CASE=h2o_co2_nc10_2d_traditional -j
make case CASE=h2o_co2_nc10_2d_benchmark -j

make run CASE=h2o_co2_nc10_2d_traditional NP=1 \
  RESULT_DIR=./results/traditional RUN_ARGS='-numSteps 1 -dt 0.01 -adaptive_dt false'
make run CASE=h2o_co2_nc10_2d_benchmark NP=1 EOS=pr \
  RESULT_DIR=./results/new-pr RUN_ARGS='-numSteps 1 -dt 0.01 -adaptive_dt false'
```

Use `EOS=sw` and `EOS=cpa` for the other new formulations. Full setup,
testing, cluster submission and recovery instructions are in
`docs/HPC_RUN.md` and `docs/CONTINUE_DEVELOPMENT.md`.
