# Common-data fit for the five-component EOS comparison

This utility fits PR, Søreide–Whitson (SW), and CPA independently to one frozen
public-data suite before the fitted packages are used in the same reservoir
case. It asks a narrow question: after each EOS has received its own best binary
parameters, which package is more consistent with the same experimental
evidence and how do the resulting flow predictions differ?

## What is held common

- components: `H2O / CO2 / CH4 / C2H6 / nC4H10`;
- raw sources: DOI-indexed NIST ThermoML records;
- data rows, preprocessing, and deterministic train/validation split;
- parameter bounds, optimizer, stopping rule, and reported metrics;
- reservoir grid, wells, initial state, schedule, and numerical tolerances.

Each EOS has its own fitted parameter vector. PR and CPA fit constant binary
interaction parameters. SW fits dry-gas binary interaction parameters and an
offset to each published aqueous correlation. Pure-component and CPA
association parameters remain fixed at their source-backed values.

The VLE records contain `T`, `P`, liquid composition `x`, and vapor composition
`y`. Their objective is the two-component log-fugacity equality residual, so a
single record simultaneously constrains the bubble branch (`P-x`), dew branch
(`P-y`), and both phase compositions. Gas-density records use the logarithmic
density ratio. The validation figure calls the resulting percentage a
**normalized residual**; it is not mislabeled as a direct bubble-pressure AARD.

## Reproduce the fit

From `tools/` under WSL or Linux:

```bash
make run-five-component-common-fit
```

The frozen input table is `data/common_calibration_data.csv`; its construction
and limitations are recorded in `data/dataset_manifest.json`. Generated fit
outputs are placed in the ignored `results/` directory.

To regenerate the bilingual figures after completing all three flow runs:

```bash
python3 example/five_component_common_fit/plot_results.py \
  --fit-results example/five_component_common_fit/results \
  --flow-root ../case/five_component_eos_tuned_compare/results/common_fit_full_final \
  --output-dir ../case/five_component_eos_tuned_compare/figures/common_fit
```

The scientific interpretation, citations, fitted values, validation results,
and flow comparison are in
`case/five_component_eos_tuned_compare/COMMON_DATA_COMPARISON.md`.
