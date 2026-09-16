# five_component_eos_tuned_legacy_compare

## Purpose

This case provides the requested traditional three-phase reference for
`five_component_eos_tuned_compare`. The comparison contains four runs:

| flow model | thermodynamic package |
|---|---|
| fully compositional O/G/W | PR |
| fully compositional O/G/W | SW |
| fully compositional O/G/W | CPA |
| Traditional oil/gas compositional + independent water | PR only |

The Traditional case has five physical components, but only
CO2/CH4/C2H6/nC4H10 enter the PR oil/gas flash. H2O has its own conservation
equation and cannot transfer into the oil or gas phase. No traditional SW or
CPA option is exposed because those models are compared only in the new fully
compositional formulation.

## Common conditions

- grid/domain: `20 x 20 x 5`, `1000 x 600 x 50 m`;
- initial state: `60 bar`, `305 K`;
- physical overall composition: H2O/CO2/CH4/C2H6/nC4H10 =
  `0.30/0.10/0.15/0.15/0.30`;
- CO2 injector: `100000 reference m3/day`, layers 0-1;
- producer: `52 bar` BHP, layers 3-4;
- schedule: 120 common targets of `0.25 day`, total `30 day`;
- closed external boundaries and identical rock, wells, TPFA transport,
  relative permeability, nonlinear tolerances and MPI rank count.

For the traditional flash, H2O is removed and the dry composition is
normalized to `CO2/CH4/C2H6/nC4H10 = 1/7, 3/14, 3/14, 3/7`. The dry PR flash
determines the oil/gas split. Its phase molar densities and the independent
water molar density are then used to obtain `Sw=0.0708563`, which reconstructs
the original global `z_H2O=0.30` instead of merely copying a saturation from a
different thermodynamic model.

## Build and run

From the repository root:

```bash
make case CASE=five_component_eos_tuned_legacy_compare -j2
make run CASE=five_component_eos_tuned_legacy_compare NP=2 \
  RESULT_DIR=./results/translated_density_30day_20260826/pr
```

The executable is deliberately PR-only; an `EOS=sw` or `EOS=cpa` selector is
neither required nor consumed.

## Reproduce figures

After completing the traditional PR run and the three fully compositional
runs:

```bash
python plot_model_comparison.py --force
python plot_field_maps.py --force
```

The script reads the volume-translated fully compositional result set at
`../five_component_eos_tuned_compare/results/translated_density_30day_20260826`
and the traditional PR output at
`results/translated_density_30day_20260826/pr`. It
writes 11 independent English and Chinese PNG/PDF figures, the plotted source
table, a phase-composition statistics table, and a provenance manifest to
`figures/translated_density_30day_20260826/`. Every model-comparison legend contains
only `Our PR`, `Our SW`, `Our CPA`, and `Traditional`.

`plot_field_maps.py` additionally writes separate vertically averaged day-30
field maps for pressure change, gas saturation, aqueous CO2 mole fraction,
gas-phase CO2 mole fraction, and water-saturation change. Each quantity uses
one shared color scale across all four models. The script does not interpolate,
smooth, normalize, or omit cells.

See [RUN_REPORT.md](RUN_REPORT.md) for conditions, numerical results,
figure-by-figure interpretation and limitations.
