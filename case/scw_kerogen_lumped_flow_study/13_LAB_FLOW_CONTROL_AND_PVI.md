# 13 — Laboratory flow-control, Darcy scaling and PVI contract

## Decision

The laboratory slab does **not** use equal injection and production rates as a physical constraint.

The control policy is:

- injector: **reservoir-volume rate control + maximum BHP**;
- producer: **fixed BHP**.

The existing Natural well-control implementation already supports this mode. A rate-controlled injector may carry a `maximumBhp` limit; after a converged state violates that limit it switches to BHP control. The producer may remain under `Bhp` control.

The old SCW kerogen benchmark used equal reservoir injection/production rates to suppress bulk pressure depletion. That is a useful legacy numerical experiment, but it is not an experimental boundary condition and is rejected for the laboratory slab.

## Pressure center

The design pressure center is

[
p_c approx 28 mathrm{MPa}.
]

This is a **reference center**, not an imposed pressure difference.

The current strict 0D PVT acceptance covers approximately 25–30 MPa. Therefore the first flow experiment must remain inside that validated window unless the 0D gate is explicitly extended and rerun.

A convenient pressure bookkeeping form is

[
p_mathrm{prod}=p_c-rac{Delta p_mathrm{design}}{2},
]

[
p_mathrm{inj,expected}=p_c+rac{Delta p_mathrm{design}}{2}.
]

This symmetric notation does not imply equal flow rates or symmetric physical boundary conditions. It merely keeps the intended pressure profile centered near 28 MPa.

The actual injector pressure is a solution of the rate-controlled flow problem.

## Rate target comes from pore-volume scale

Let

- (PV_mathrm{eff}) be the measured accessible pore volume at operating T/P;
- (r_mathrm{PVI}) be the desired injected pore volumes per second.

Then the injector reservoir-volume rate target is

[
Q_mathrm{inj,target}
=
r_mathrm{PVI}PV_mathrm{eff}.
]

The corresponding nominal residence time is

[
t_mathrm{res}
=
rac{PV_mathrm{eff}}{Q_mathrm{inj,target}}
=
rac{1}{r_mathrm{PVI}}.
]

Thus the experimental rate should be selected by specifying a meaningful residence/PVI scale, not by copying the old `0.25 PV/day` numerical benchmark.

## Darcy velocity

For measured effective open flow area (A),

[
u_D=rac{Q}{A}.
]

For measured effective porosity (phi_mathrm{eff}), an interstitial/pore-velocity estimate is

[
u_papproxrac{Q}{phi_mathrm{eff}A}.
]

Both quantities should be reported because a pump flow rate alone is not transferable between slab geometries.

## First pressure-drop estimate

Before any multiphase simulation is used to set well pressures, perform a single-phase Darcy estimate:

[
Delta p_mathrm{sp}
=
rac{mu_mathrm{ref}LQ}
     {k_mathrm{abs}A}.
]

All quantities must correspond to the actual experiment:

- (L): measured effective flow length;
- (A): measured effective open flow area;
- (k_mathrm{abs}): measured specimen permeability, preferably checked at operating T/P;
- (mu_mathrm{ref}): validated target-T/P viscosity for the intended initial/injected single-phase reference;
- (Q): proposed experimental rate.

This is an order-of-magnitude design relation, not a fitted multiphase pressure drop.

Only after this estimate is available should (p_mathrm{prod}), expected injector pressure and BHP cap be fixed.

## Multiphase correction

After the single-phase scale is understood, a mobility estimate may be formed:

[
lambda_t=sum_alpharac{k_{r,alpha}}{mu_alpha},
]

[
Delta p_mathrm{mp}
sim
rac{QL}{k_mathrm{abs}Alambda_t}.
]

Because the first-stage relative-permeability curves are Corey sensitivity assumptions rather than SCW calibrations, this multiphase estimate is also a sensitivity range.

It must not be used to retroactively tune the Corey curves to a desired pressure drop.

## Injector maximum BHP

The maximum injector BHP is a constraint, not the nominal driving pressure.

It must satisfy all relevant limits:

[
p_mathrm{inj,max}
=
min(
p_mathrm{apparatus,safe},
p_mathrm{PVT,validated,max},
p_mathrm{experiment,design cap}
).
]

No numerical value is committed until the apparatus pressure rating, pump limit and Darcy-scale calculation are documented.

With the current 0D gate, exceeding 30 MPa also requires an extended zero-dimensional PVT validation.

## Producer BHP

The producer is fixed-BHP.

Its target is selected after the Darcy calculation so that:

- the pressure drop is measurable but not excessive;
- the entire expected pressure field remains in the accepted PVT range;
- the pump/control hardware can maintain the requested inlet rate without immediately hitting maximum BHP;
- residence time remains in the intended experimental range.

The producer BHP is therefore not chosen merely as “initial pressure minus 0.5 bar.”

## No equal-volume injection/production constraint

Under the laboratory policy:

- injection volume is prescribed while rate control is active;
- production rate is an output of fixed outlet BHP, fluid compressibility, phase behavior, storage and transport;
- total in-place inventory may temporarily change;
- if the injector reaches its maximum BHP, its actual rate changes.

This is the experimentally meaningful behavior.

The mass-balance/inventory diagnostics remain mandatory.

## PVI is the primary time coordinate

All final flow comparisons should use cumulative **actual** injected pore volumes:

[
mathrm{PVI}(t)
=
rac{
int_0^t Q_{mathrm{inj,actual,res}}(	au),d	au
}
{PV_mathrm{eff}}.
]

The word `actual` is important. If the injector switches from rate control to maximum-BHP control, nominal target rate no longer represents injected volume.

Therefore plots/tables should use:

- primary x-axis: PVI;
- secondary/reporting coordinate: physical time in s/min/h as appropriate.

“Days” must not be the sole normalization for an experimental slab.

## Paired 360/380 °C control rule

For the formal SCW control experiment, the target dimensionless injection intensity is identical:

`r_PVI,360 = r_PVI,380`.

The physical rate is derived independently from the measured operating-condition pore volume:

`Q(T) = r_PVI * PV_eff(T,p)`.

The pressure drop is **not** forced to be equal between the two temperatures. The difference in injector BHP and slab pressure drop at matched PVI is an experimental response because density, viscosity, phase split and mobility are temperature dependent.

Likewise, the producer BHP/control rule and injector maximum-BHP rule are identical between the pair; they are not retuned separately to force equal production.

## Required outputs

At minimum record:

- physical time;
- cumulative actual injected reservoir volume;
- PVI;
- injector active control mode;
- injector actual reservoir rate;
- injector BHP;
- producer BHP;
- producer actual reservoir rate;
- total pressure drop;
- produced component masses;
- domain component inventory and closure.

This makes control switching and storage effects auditable.

## Current status

The control **policy** is accepted, but the numerical targets are intentionally not filled.

Missing inputs are:

- actual slab (L,A);
- (PV_mathrm{eff}(T,p));
- (k_mathrm{abs}(T,p));
- validated reference viscosity;
- pump-rate range;
- apparatus safe pressure;
- selected target residence/PVI rate.

Until those exist:

`FLOW_CONTROL_GATE_BLOCKED`.

This is a design-data block, not a numerical-solver block.
