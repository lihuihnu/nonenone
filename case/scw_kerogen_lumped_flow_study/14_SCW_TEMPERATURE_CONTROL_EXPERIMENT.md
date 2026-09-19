# 14 — Formal 360 °C / 380 °C SCW temperature-control experiment

## Primary scientific question

The formal flow experiment must answer a paired question:

> What changes when the water-containing system crosses the critical temperature of water while pressure, composition, porous-medium geometry and flow-control design are otherwise held fixed?

The mandatory pair is:

- **control:** 360 °C / 28 MPa;
- **SCW test:** 380 °C / 28 MPa.

Water's critical point is approximately 373.946 °C and 22.064 MPa. At 360 °C the temperature is below the critical temperature; because 28 MPa is above the saturation pressure for a subcritical state this close to the critical point, the water reference state is compressed liquid. At 380 °C and 28 MPa, both temperature and pressure exceed the critical values, so the water reference state is supercritical.

The experiment is therefore a controlled crossing of the water critical temperature at the same nominal pressure.

## Near-critical diagnostic

A third state is retained:

- 374 °C / 28 MPa.

This is a **diagnostic**, not the primary control.

Its purpose is to show whether the 360 -> 380 response is approximately smooth or becomes strongly nonlinear in the near-critical region.

The core causal comparison remains 360 vs 380.

## What must remain identical

The formal pair shares:

- the same 60×20×1 topology;
- the same apparatus dimensions and inlet/outlet manifolds;
- the same porous-medium batch/specification;
- the same initial overall composition;
- the same injection composition;
- the same 28 MPa pressure center;
- the same target PVI rate;
- the same injector rate + maximum-BHP control logic;
- the same fixed-BHP producer logic;
- the same gravity orientation;
- the same output schedule in PVI;
- the same EOS identity within each PR or CPA paired comparison.

The only prescribed state-variable change is temperature.

The machine-readable definition is \`porous_media/scw_temperature_control_pair.csv\`.

## Two levels of interpretation

### Level 1 — mechanism-isolation numerical pair

For the first numerical comparison, the rock-fluid closures are deliberately frozen between 360 and 380 °C:

- the same Corey sensitivity curve;
- the same residual saturations;
- the same capillary-pressure sensitivity case;
- the same homogeneous porosity/permeability field.

This isolates the consequences of temperature acting through:

- phase equilibrium;
- phase composition;
- density;
- viscosity;
- phase mobility through the viscosity term;
- component partitioning and transport.

No Corey exponent, residual saturation or capillary parameter may be retuned separately at 380 °C to improve the appearance of the SCW case.

### Level 2 — full laboratory temperature response

A real porous medium may itself respond to temperature.

If later measurements show temperature-dependent:

- absolute permeability;
- wettability;
- interfacial tension;
- capillary pressure;
- relative permeability;
- accessible pore volume,

those data belong to a separately labelled **full-physics paired analysis**.

They are part of the total experimental temperature effect, but they must not be silently inserted into the Level-1 mechanism-isolation comparison.

This separation lets the project distinguish:

1. fluid thermodynamic/transport response to crossing the critical state;
2. porous-medium/rock-fluid response to the temperature change.

## Flow-rate matching

The formal pair uses the same target **PVI rate**:

\[
r_{\mathrm{PVI},360}
=
r_{\mathrm{PVI},380}.
\]

For each run,

\[
Q_\mathrm{inj,target}(T)
=
r_\mathrm{PVI}PV_\mathrm{eff}(T,p).
\]

If measured \(PV_\mathrm{eff}\) changes slightly with temperature, the physical volumetric rate may therefore differ correspondingly. This is not an independent tuning variable; it is the consequence of holding the dimensionless residence/PVI scale fixed.

A supplementary same-\(Q\) engineering comparison may be reported, but it does not replace the primary matched-PVI experiment.

## Pressure-control matching

Both runs use the same design pressure center:

\[
p_c=28\ \mathrm{MPa}.
\]

The same control law is used in both experiments:

- injector: reservoir-volume rate control with maximum BHP;
- producer: fixed BHP.

The target pressure drop is not separately tuned at each temperature.

Because viscosity changes strongly with temperature, the injector BHP and slab pressure drop are expected outcomes and are important observables.

## 0D prerequisite at 28 MPa

The formal pair requires a dedicated 0D anchor before flow:

\[
z=z_\mathrm{BASE},\qquad
z_\mathrm{H2O}=0.20,\qquad
p=28\ \mathrm{MPa}
\]

at both:

\[
T=633.15\ \mathrm{K}
\]

and

\[
T=653.15\ \mathrm{K}.
\]

PR and CPA must each independently pass:

- flash convergence;
- final stability;
- phase-role consistency;
- material closure;
- fugacity closure;
- finite positive density and viscosity.

The existing broad 25–30 MPa scan contains these states, but the formal experiment now records them as named paired anchors rather than relying only on their membership in a large scan.

## Current 0D paired-anchor result

The dedicated 28 MPa paired anchor has now run successfully for both EOS:

\`SCW_360_380_CONTROL_PAIR_28MPA = PASS\`.

At the current BASE composition \(z_{\mathrm{H2O}}=0.20\), all three 28 MPa anchor temperatures remain a single Oil-role phase in both PR and CPA.

The current 0D values are:

| EOS | T °C | phase count | mass density kg/m³ | LBC viscosity Pa·s |
|---|---:|---:|---:|---:|
| PR | 360 | 1 | 681.648 | 1.97572e-4 |
| PR | 374 | 1 | 676.596 | 1.92852e-4 |
| PR | 380 | 1 | 674.406 | 1.90867e-4 |
| CPA | 360 | 1 | 614.440 | 1.45707e-4 |
| CPA | 374 | 1 | 610.172 | 1.43371e-4 |
| CPA | 380 | 1 | 608.323 | 1.42380e-4 |

For the formal 360 -> 380 comparison at the initial BASE state:

- PR density changes by approximately \(-1.06\%\);
- PR viscosity changes by approximately \(-3.39\%\);
- CPA density changes by approximately \(-1.00\%\);
- CPA viscosity changes by approximately \(-2.28\%\).

This initial-state result does **not** show a phase-topology change.

That is scientifically important: the formal experiment must not equate “water crosses its critical point” with “the hydrocarbon-rich mixture must immediately split into a new phase.” The stronger SCW signal may instead emerge as injected water drives the local overall composition into water-rich parts of the phase diagram.

Therefore the flow study must compare matched-PVI trajectories and water-rich local states, not only the initial bulk state.

## Primary observables

The required observables are registered in \`porous_media/scw_temperature_control_observables.csv\`.

### Thermodynamic observables

At matched pressure/composition:

- phase count and phase role;
- phase fractions;
- phase compositions;
- density;
- viscosity.

### Transport observables

At matched PVI:

- phase mobility;
- slab pressure drop;
- injector BHP;
- phase saturation fields;
- front position;
- breakthrough PVI.

### Component-selective transport

For each conserved pseudo-component \(i\), define

\[
RF_i(\mathrm{PVI})
=
\frac{M_{i,\mathrm{produced,cumulative}}}
     {M_{i,\mathrm{initial\ in\ place}}}.
\]

A simple effluent enrichment metric is

\[
E_i
=
\frac{w_{i,\mathrm{produced}}}
     {w_{i,\mathrm{initial}}},
\]

reported only where the denominator is nonzero.

Pairwise selective transport may be summarized as

\[
S_{ij}
=
\frac{RF_i}{RF_j}.
\]

The temperature effect is then compared at the same PVI, for example

\[
\Delta RF_i
=
RF_{i,380}-RF_{i,360}.
\]

This directly addresses whether light, middle and heavy pseudo-components are transported differently after crossing the water critical state.

Any breakthrough definition must preregister its threshold before looking at the final curves.

## Producer composition is a primary paired endpoint

The formal 360/380 comparison must consume the first-class producer output defined in \`15_PRODUCER_COMPOSITION_AND_SELECTIVITY.md\`.

At every matched PVI the analysis separates:

1. total hydrocarbon recovery;
2. component-resolved recovery fractions \(RF_i\);
3. instantaneous producer mass fractions \(Y_i\);
4. instantaneous and cumulative \(E_{L/H}\).

A higher total recovery at 380 °C is not, by itself, evidence for selective SCW transport. Selectivity requires a component-resolved change after normalization to the initial in-place composition.

## Experimental pairing and history control

A formal laboratory claim requires control of specimen history.

Preferred practice is:

- specimens cut/packed from the same characterized batch;
- matched \(\phi\), \(k\) and PV before the run;
- identical feed preparation;
- identical pressure and flow-control protocol;
- record run order and thermal exposure;
- do not reuse a specimen after 380 °C SCW exposure as the nominally pristine 360 °C control unless post-exposure characterization proves no material change.

Replicate strategy and run order are apparatus/experimental-design inputs and remain to be preregistered.

## EOS comparison rule

PR and CPA are each evaluated with the same 360/380 pair.

The scientific comparisons are therefore:

\[
\mathrm{PR}_{380}-\mathrm{PR}_{360}
\]

and

\[
\mathrm{CPA}_{380}-\mathrm{CPA}_{360}.
\]

Only after those within-EOS temperature effects are established should the project compare whether PR and CPA predict different magnitudes or mechanisms.

A visually appealing 380 °C flow result is not evidence that an EOS is more physical.

## Gate

The formal temperature-control gate is in:

\`porous_media/scw_temperature_control_gate.csv\`.

The pair definition is accepted, but flow remains blocked until:

- the paired 28 MPa 0D anchors are explicitly recorded;
- property validation is admissible at both temperatures;
- the actual slab/apparatus and flow-control gates pass;
- paired output metrics and mass balance are implemented.
