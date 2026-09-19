# 11 — Laboratory-reproducible porous-medium baseline

## Decision

The next flow model must not be described as an abstract reservoir.

The porous medium is to be defined as a **laboratory specimen** with measurable dimensions, porosity, absolute permeability, wettability, pore volume and multiphase constitutive behavior.

For the first experimental stage, the preferred reference is **fired Berea sandstone**. A 2-D packed-sand slab remains an attractive second route for visualization and repeatability, but every packing must be treated as a new specimen whose porosity, permeability and effective pore volume are measured rather than copied from literature.

The previous `scw_kerogen_common` values

- `phi = 0.35`;
- `kx = ky = 1500 mD`;
- `kz = 150 mD`;
- `L = 1.20 m`;

are now classified as **legacy numerical values only**. They are not admissible as the porous-medium definition of the laboratory-reproducible SCW study.

## Why fired Berea is the first reference

Niu, Al-Menhali and Krevor used a homogeneous fired Berea core with:

- diameter = 3.8 cm;
- length = 20 cm;
- measured porosity = 21%;
- water permeability = 212 mD;
- pore volume = 47.6 mL.

The core was heated at 700 °C to stabilize mobile clays before repeated multiphase experiments. Berea is also the basis of classic water-wet three-phase relative-permeability studies.

This gives a much stronger experimental anchor than assigning an arbitrary sandstone permeability to a numerical grid.

These values are **reference-specimen values**, not values that may be copied into a future specimen without measurement.

## Candidate porous media

See `porous_media/media_candidates.csv`.

### Fired Berea sandstone

Preferred for the first coreflood because it has:

- standard laboratory supply;
- moderate permeability, which is advantageous when SCW viscosity is very low;
- extensive water-wet two- and three-phase flow literature;
- capillary-pressure data;
- a published clay-stabilization procedure.

### Bentheimer sandstone

Bentheimer is highly homogeneous and quartz-rich, with typical porosity around 0.21–0.27 and permeability roughly 0.5–3 D.

It is a good reproducibility benchmark, but the high permeability can make pressure-drop measurements difficult for a very low-viscosity SCW phase.

### Artificial water-wet sand pack

Kianinejad et al. used 0.25 mm sand with `phi=0.365` and `k=13 D`; Moghadasi et al. report another water-wet sand-pack reference with `phi≈0.37` and `k≈2.9 D`.

These are valuable packing references, but neither permeability should be copied into a future slab. The slab permeability is an **experimental output of the packing protocol**.

If a 2-D slab is selected, grain-size distribution, packing vibration/compaction, confining stress and thermal cycle become part of the reproducibility specification.

## Porosity, permeability and effective pore volume

The simulation must consume the properties of the actual specimen.

Minimum measurements are listed in `porous_media/measurement_plan.csv`.

The acceptance sequence is:

1. measure physical dimensions and bulk volume;
2. measure dry mass / mineral or grain identity;
3. measure ambient effective porosity and pore volume;
4. measure ambient absolute permeability by single-phase flow;
5. repeat single-phase permeability under the target effective stress and, where feasible, at 360/374/380 °C and 25–30 MPa;
6. determine accessible pore volume under operating conditions;
7. repeat key measurements after SCW exposure.

The effective pore volume used for injection rates must therefore be

`PV_effective(T,p) = measured accessible pore volume`,

not merely `phi_literature * numerical bulk volume`.

The published Berea `47.6 mL` is a useful scale for a 38 mm × 200 mm core, but it is not a universal pore volume.

## Wettability

The current prior is **water-wet**, because fired Berea and clean quartz-rich sandstones are generally water-wet in the lower-temperature analogue literature.

This is not a 380 °C SCW calibration.

At the target temperatures the oil/water/rock surface chemistry, mineral dissolution and interfacial tension can differ markedly. Therefore wettability is a hard experimental/sensitivity item.

Until target-condition evidence exists:

- water-wet is the reference prior;
- neutral/mixed-wet behavior must be included as a sensitivity interpretation;
- no numerical contact angle is labelled measured.

## Relative permeability

There is no adequate direct relative-permeability dataset for the present SCW + experiment-driven oil pseudo-component system at 360–380 °C and 25–30 MPa.

Therefore the first flow stage may use only a **simplified Corey sensitivity family**.

The committed scenarios are in `porous_media/corey_sensitivity.csv`.

The reference-proxy case combines lower-temperature Berea information:

- `Swr = 0.1699`;
- `Sor = 0.18`;
- `Sgr = 0` for drainage;
- `nw = 6`;
- `no = 2`;
- `ng = 3`;
- normalized Corey endpoints = 1.

This is deliberately labelled `SENSITIVITY_ONLY`.

It is **not one experimental dataset**: water/gas quantities come from a CO2/brine Berea design and oil residual/exponent information comes from historical three-phase/oil-water Berea studies.

Low- and high-trapping scenarios must accompany it so that any SCW-flow conclusion is not conditional on one arbitrary Corey curve.

If target-condition relative permeability is later measured, these screening curves are replaced rather than tuned to match the desired flow result.

## Residual saturations

Literature demonstrates substantial spread even within Berea.

Examples used only to define the sensitivity envelope include:

- `Swr ≈ 0.11` in one Berea CO2/brine study;
- `Swr = 0.1699` in the Niu experimental-design model;
- historical Berea connate-water saturations around 0.32–0.36;
- Oak-derived three-phase residual-oil fits around 0.08 and 0.18;
- historical waterflood Berea residual oil around 0.24–0.31.

This variation is exactly why residual saturation is treated as a sensitivity parameter, not a universal Berea constant.

## Capillary pressure

Capillary pressure is potentially important at laboratory scale.

Low-temperature Berea references include Brooks-Corey proxies such as:

- `Pe = 2.5 kPa, lambda = 0.67`;
- `Pe = 5.5 kPa, lambda = 0.6`.

MICP can give a much larger mercury-air entry pressure (around 53.5 kPa in another Berea study), illustrating that **entry pressure is fluid-pair dependent**.

Therefore:

- no low-temperature `Pe` is called an SCW-hydrocarbon calibration;
- `Pc=0` is retained only as a control calculation;
- Berea proxy curves are sensitivity cases;
- final use should rely on target-condition IFT/wettability or a directly measured displacement curve.

A critical implementation finding is that the current Natural face-flux path uses a common pressure for all phases; an explicit phase capillary-pressure closure is not presently part of the flow kernel.

Consequently the porous-medium gate cannot pass until either:

1. capillary pressure is implemented and tested, or
2. a documented force-scale analysis demonstrates that `Pc=0` is an acceptable experimental regime for the selected flow rate and specimen.

Because SCW has low viscosity and published Berea entry pressures are on the kPa scale, this issue should not be dismissed automatically.

## Pseudo-3D slab discretization

The laboratory slab uses a fixed `60×20×1` regular Cartesian topology while leaving all physical dimensions apparatus-defined.

See `12_LAB_PSEUDO3D_SLAB_DESIGN.md` and `porous_media/pseudo3d_slab_geometry_contract.csv`.

Only the topology is fixed. The old numerical benchmark dimensions and point-well locations are not inherited.

Stage 1 is spatially homogeneous by design; stochastic permeability/porosity fields are deferred until the homogeneous experiment has been validated.

Because `nz=1`, the slab is thickness-averaged. It can resolve in-plane buoyancy if a model in-plane axis is aligned with physical gravity, but it cannot resolve or support claims about true out-of-plane/3-D vertical gravity segregation.

## 2-D slab strategy

For the eventual visual 2-D experiment, two defensible routes exist.

### Route A — sandstone slab

Cut a rectangular slab from the same characterized Berea block. Measure both in-plane permeability directions.

Advantages:

- mineralogy and pore structure are physically sandstone;
- direct link to the coreflood reference.

Disadvantages:

- custom high-P/T sealing;
- slab-to-slab heterogeneity;
- anisotropy must be measured.

### Route B — packed quartz-sand slab

Use a specified sieve distribution and reproducible compaction protocol inside a rectangular high-P/T cell.

Advantages:

- easy geometry replication;
- pack can be rebuilt;
- useful for direct visualization.

Disadvantages:

- porosity/permeability depend on each packing;
- high permeability may lead to extremely small viscous pressure drop;
- grain rearrangement and SCW mineral-surface changes must be checked.

The numerical model must inherit the **measured slab properties**, not the properties of a literature sand pack.

## Flow-entry gate

The porous-medium gate is in `porous_media/porous_media_entry_gate.csv`.

The existing zero-dimensional PR/CPA PVT gate has passed, but this new porous-medium gate is blocked by design until an actual laboratory specimen/cell is selected and characterized.

The global flow comparison therefore remains frozen.
