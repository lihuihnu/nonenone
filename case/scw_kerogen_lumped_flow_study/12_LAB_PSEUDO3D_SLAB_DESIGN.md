# 12 — Experimental-scale 60×20×1 pseudo-3D slab design

## Design decision

The `60×20×1` mesh is retained, but it is no longer assigned reservoir-scale dimensions.

It is a **fixed numerical topology for a laboratory slab**:

[
N_x=60,qquad N_y=20,qquad N_z=1.
]

The physical dimensions are determined by the actual high-temperature/high-pressure experimental cell:

[
L_x=L_{x,mathrm{apparatus}},quad
L_y=L_{y,mathrm{apparatus}},quad
L_z=L_{z,mathrm{apparatus}}.
]

Only after those dimensions are measured are the uniform cell dimensions defined:

[
Delta x=rac{L_x}{60},qquad
Delta y=rac{L_y}{20},qquad
Delta z=L_z.
]

This prevents the numerical model from dictating the apparatus size.

The machine-readable contract is `porous_media/pseudo3d_slab_geometry_contract.csv`.

## What “pseudo-3D” means here

The slab has a real physical thickness (L_z), but the numerical model contains only one cell through that thickness.

Therefore each cell represents a thickness-averaged control volume:

[
V_mathrm{cell}
=
rac{L_xL_yL_z}{60	imes20}.
]

The model can still carry three Cartesian permeability/transmissibility directions in its numerical machinery, but it does **not** spatially resolve gradients through the thickness.

The term `pseudo-3D slab` is used only to indicate a finite-thickness laboratory body represented by a 2-D in-plane grid plus one thickness-averaged layer.

It must not be described as a fully resolved 3-D porous-medium experiment.

## Apparatus-first geometry

No final values of (L_x,L_y,L_z) are committed yet.

The future apparatus/cell design determines:

- usable flow length;
- active porous width or height;
- packed/rock thickness;
- inlet and outlet manifold geometry;
- observation-window geometry if used;
- pressure-tap positions.

The numerical model then inherits those dimensions.

A convenient integer cell size may be considered when choosing an apparatus, but it is not a scientific requirement and must not override pressure rating, sealing, thermal control, optical access or specimen constraints.

For example, if a future apparatus happens to have dimensions that divide cleanly by 60 and 20, the resulting (Delta x,Delta y) will be convenient. That is a meshing benefit, not a basis for choosing the apparatus dimensions.

## Stage-1 homogeneous porous medium

The first flow experiment is deliberately homogeneous.

For every cell:

[
phi_{ijk}=phi_mathrm{measured}
]

and, for an isotropic first approximation,

[
k_{x,ijk}=k_{y,ijk}=k_mathrm{measured}.
]

If the slab material has measured in-plane anisotropy, a constant tensor may be used:

[
k_x
eq k_y,
]

but the value remains spatially uniform.

No random permeability multiplier, geostatistical realization, layered facies, channel or porosity noise is allowed in Stage 1.

The reason is experimental identifiability: the first comparison is intended to isolate

- PR vs CPA thermodynamics;
- phase appearance/disappearance;
- density and viscosity;
- relative-permeability sensitivity;
- capillary-pressure sensitivity;
- injection-rate and gravity effects.

Adding geological heterogeneity before these effects are validated would make model discrepancies non-identifiable.

The staged policy is stored in `porous_media/heterogeneity_stage_policy.csv`.

## Future heterogeneity

Random or image-derived permeability fields are deferred.

They become scientifically meaningful only after the homogeneous slab has a validated experimental baseline.

A future stochastic field must specify at least:

- distribution form;
- arithmetic/geometric mean;
- variance;
- spatial correlation lengths;
- anisotropy;
- random seed or image provenance;
- whether porosity is correlated with permeability.

The field must preserve the experimentally measured bulk (phi), effective permeability and accessible pore volume within documented tolerances.

A random field is therefore a later hypothesis test, not a first-stage realism decoration.

## Gravity and the nz=1 limitation

The gravity scope is explicitly limited.

### Gravity-off reference

A gravity-off case is valid and useful:

[
mathbf g=mathbf 0.
]

It isolates pressure-driven transport and EOS/property effects.

### Vertical slab with gravity in the model plane

If the physical slab is mounted vertically and the model (y)-axis is defined as vertical, gravity may be applied in-plane.

Then the (60	imes20) grid can resolve a **2-D in-plane buoyancy segregation pattern** over the 20 cells in the gravity direction.

This is a valid laboratory sensitivity study.

However, it is still not a fully three-dimensional gravity-segregation experiment.

### Out-of-plane gravity

If physical vertical is normal to the slab plane, the model z direction has only one cell:

[
N_z=1.
]

No saturation or composition gradient can be resolved through the slab thickness.

Therefore this model cannot be used to claim:

- resolved vertical gravity override through thickness;
- vertical fingering through multiple layers;
- 3-D buoyant plume rise;
- genuine vertical gravitational segregation in the unresolved z direction.

Those questions require:

[
N_z>1
]

and a physical apparatus with a resolved vertical dimension.

The allowed gravity claims are machine-readable in `porous_media/gravity_scope.csv`.

## Grid-quality rule

The grid is regular Cartesian and uniform:

[
Delta x = L_x/60,qquad
Delta y = L_y/20,qquad
Delta z=L_z.
]

This gives exactly 1200 cells.

The first-stage grid must not be locally refined around inlet/outlet wells merely to improve numerical appearance. If the apparatus uses distributed manifolds rather than point ports, the boundary condition should reproduce that geometry instead of hiding it through grid refinement.

A grid-convergence check is still required later. A suitable comparison is, for the same apparatus and physics:

- `30×10×1`;
- `60×20×1`;
- `120×40×1`.

The production reference remains `60×20×1` only if integral recovery, pressure drop and front position are sufficiently grid-independent.

## Effective pore volume and injection-rate scaling

The numerical bulk volume is

[
V_b=L_xL_yL_z.
]

But injection rates in pore volumes must use the **measured operating-condition effective pore volume**:

[
PV_mathrm{eff}=PV_mathrm{measured}(T,p),
]

not automatically

[
phi V_b.
]

The latter is a consistency check only.

Thus a target injection rate expressed as pore volumes per unit time becomes

[
Q = mathrm{PVI rate}	imes PV_mathrm{eff}.
]

This keeps the numerical rate tied to the experimental specimen rather than to an assumed porosity.

## Relationship to the existing 60×20×1 code

The legacy `scw_kerogen_common/benchmark_2d_common.hpp` also uses `60×20×1`, but its dimensions and rock values are historical numerical choices:

- (L_x=1.2) m;
- (L_y=0.1) m;
- (L_z=0.1) m;
- (phi=0.35);
- (k_x=k_y=1500) mD;
- (k_z=150) mD.

Those values must not be inherited by the new laboratory slab.

Only the **cell-count topology** `60×20×1` is retained.

A future executable laboratory case should read or compile in apparatus/specimen measurements explicitly and should fail if those values remain `TBD`.

## Well-control coupling

The apparatus geometry does not merely set cell dimensions; it also sets the experimental flow scale.

Once (L_x,L_y,L_z), effective open area, measured (PV_{m eff}), (k) and validated viscosity are available, the well-control design follows `13_LAB_FLOW_CONTROL_AND_PVI.md`.

The slab does not use equal injection and production rates.

- injector: rate control with maximum-BHP limit;
- producer: fixed BHP;
- target rate: derived from PVI/residence-time scale;
- pressure drop: checked first with a Darcy order-of-magnitude estimate;
- final comparison coordinate: actual cumulative PVI.

## Entry-gate effect

The geometry concept itself is accepted:

`pseudo3d_grid_topology = PASS_DESIGN`.

But the physical geometry remains blocked:

`apparatus_dimensions = OPEN`.

Therefore no new flow case is allowed until the actual slab/core holder provides (L_x,L_y,L_z), the corresponding porous specimen has measured (phi,k,PV_mathrm{eff}), and the other porous-medium gates are satisfied.
