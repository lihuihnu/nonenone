# 15 — Producer composition as a first-class SCW output

## Decision

Producer composition is a **primary scientific output**, not a post-processing convenience.

Every formal 360/380 °C run must write a dedicated:

\`producer_composition.csv\`

from the exact conserved-component well source.

The output must not reconstruct composition from saturation, phase labels or a sampled cell next to the producer.

## Primary quantities

For every conserved component \(i\), the producer output records:

\[
\dot m_{i,\mathrm{prod}}
=
\max(0,-\dot m_{i,\mathrm{well}})
\]

and the instantaneous produced mass fraction

\[
Y_i
=
\frac{\dot m_{i,\mathrm{prod}}}
{\sum_j \dot m_{j,\mathrm{prod}}}.
\]

The cumulative component production is integrated at every accepted internal time step:

\[
M_{i,\mathrm{prod}}(t)
=
\sum_n
\dot m_{i,\mathrm{prod}}^{n+1}\Delta t_n.
\]

This uses the same accepted end-step Backward-Euler convention as the component mass-balance ledger.

## Recovery fraction

For every hydrocarbon lump:

\[
RF_i(t)
=
\frac{M_{i,\mathrm{prod}}(t)}
{M_{i,0}},
\]

where \(M_{i,0}\) is the initial global in-place conserved-component mass.

H2O may also be written with the same mathematical quantity for completeness, but hydrocarbon recovery fractions are the primary recovery metrics because injected H2O can exceed the initial in-place H2O inventory.

Total hydrocarbon recovery must be reported separately from the individual-lump recovery fractions.

## Light/heavy enrichment

The enrichment diagnostic uses a consistent **mass basis**.

For the current experiment-driven four-oil-lump topology, define the derived Light group as:

\[
L =
\mathrm{OIL\_GASOLINE}
+
\mathrm{OIL\_DIESEL},
\]

and:

\[
H =
\mathrm{OIL\_HEAVY}.
\]

Middle remains an independent primary output and is excluded from the L/H diagnostic.

The instantaneous lightening indicator is:

\[
E_{L/H}
=
\frac{Y_L/Y_H}
     {w_{L,0}/w_{H,0}}.
\]

This is intentionally not written as a mass-fraction numerator divided by a mole-fraction \(z_L/z_H\) denominator. If \(z\) denotes mole fractions, mixing those bases would create a nonphysical enrichment number.

The cumulative analogue is:

\[
E_{L/H,\mathrm{cum}}
=
\frac{M_{L,\mathrm{prod}}/M_{H,\mathrm{prod}}}
     {M_{L,0}/M_{H,0}}.
\]

Interpretation:

- \(E_{L/H}>1\): producer stream is enriched in the registered light group relative to the initial in-place light/heavy mass ratio;
- \(E_{L/H}=1\): no net L/H enrichment on that basis;
- \(E_{L/H}<1\): relative heavy enrichment.

No threshold is used as a model score. The full trajectory versus PVI is retained.

## Why recovery and selectivity are separate

Two physically different mechanisms must not be conflated.

### Higher total recovery

SCW may increase total produced hydrocarbon mass without changing the relative composition much.

That is a recovery/mobility effect.

### Selective light-component transport

SCW may preferentially mobilize or partition lower-boiling components so that:

\[
RF_L-RF_H
\]

or

\[
E_{L/H}
\]

changes even when total hydrocarbon recovery changes little.

That is a component-selectivity effect.

The formal 360/380 analysis therefore reports both:

1. total hydrocarbon recovery;
2. component-resolved \(RF_i\);
3. instantaneous and cumulative \(E_{L/H}\).

## PVI alignment

Each producer-composition row also stores cumulative actual injected reservoir volume and:

\[
\mathrm{PVI}(t)
=
\frac{\int_0^t Q_{\mathrm{inj,actual,res}}\,dt}
{PV_\mathrm{eff}}.
\]

For the formal laboratory slab, \(PV_\mathrm{eff}\) must be the measured operating-condition effective pore volume.

If no measured PV has been configured, the output leaves PVI undefined rather than substituting a geometric or literature pore volume.

## Component-role boundary

The current experiment-driven 5-component PVT fluid is:

- H2O;
- OIL_GASOLINE;
- OIL_DIESEL;
- OIL_MIDDLE;
- OIL_HEAVY.

It does **not** contain an independent gas-product lump.

Therefore the project must not rename:

- the Gas phase as a conserved Gas component; or
- OIL_GASOLINE as Gas.

If a future formal flow fluid requires the requested explicit output set:

- H2O;
- Gas;
- Light;
- Middle;
- Heavy,

the Gas component must first be added through the already-required unified gas + recovered-liquid mass/mole basis, with no double-counting of low-boiling losses.

Until then, the first-class output uses the actual conserved component names and writes the derived Light group only as a selectivity diagnostic.

## Implementation

The generic producer ledger is:

\`output/include/output/metrics/producer_composition.hpp\`.

The first-class CSV writer is:

\`output/include/output/well/producer_composition_output.hpp\`.

The common case runner integrates it at every accepted internal step.

The output contains, for every producer and every fixed output state:

- physical time;
- cumulative actual injected reservoir volume;
- PVI;
- instantaneous producer component mass rate;
- instantaneous component mass fraction;
- cumulative produced component mass;
- component recovery fraction;
- instantaneous mass-basis \(E_{L/H}\);
- cumulative mass-basis \(E_{L/H}\).

The exact field contract is machine-readable in:

\`porous_media/producer_composition_output_contract.csv\`.

## Runtime smoke evidence

The first-class output path has been compiled and executed through a real SCW kerogen reservoir case.

A one-accepted-step \`scw_kerogen_lmh_1d\` smoke run generated \`producer_composition.csv\` from the production well-source kernel.

The accepted row satisfied:

- instantaneous producer mass fractions sum to 1;
- cumulative component production is finite and non-negative;
- each \(RF_i\) is finite for components with positive initial inventory;
- \`RF_total_hydrocarbon\` is finite and non-negative;
- instantaneous and cumulative \(E_{L/H}\) are finite and positive.

For the homogeneous short-step regression, all hydrocarbon component recovery fractions and total hydrocarbon recovery are approximately \(2.50006\times10^{-5}\), while \(E_{L/H}=1\). This is the expected no-selectivity control behavior and demonstrates that the metric does not manufacture a lightening signal when components are produced proportionally.

The formal laboratory PVI remains undefined in this legacy smoke because no measured operating-condition \(PV_\mathrm{eff}\) is supplied. The writer intentionally leaves PVI as NaN rather than inserting a geometric or literature pore volume.

## Formal paired comparison

At matched PVI, the 360/380 comparison must report:

\[
\Delta RF_i
=
RF_{i,380}-RF_{i,360},
\]

and:

\[
R_{E,L/H}
=
\frac{E_{L/H,380}}
     {E_{L/H,360}},
\]

where both are defined.

This lets the study distinguish a change in total displacement efficiency from a change in compositional selectivity.

The formal gate is:

\`porous_media/producer_composition_gate.csv\`.
