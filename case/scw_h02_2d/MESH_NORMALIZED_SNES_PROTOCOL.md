# H02 mesh-normalized SNES convergence protocol

Registered before changing the H02 nonlinear convergence implementation.

## Trigger condition

This protocol is activated only if the exact 60-s-aligned 60x20 B failure
confirms the pattern already seen in the auxiliary short run:

- global scaled residual L2 norm remains slightly above 1e-6;
- residual RMS is O(1e-8);
- residual Linf is well below 1e-6;
- no single conservation/equilibrium/closure/well row is itself O(1e-6).

That pattern means the fixed global L2 absolute tolerance becomes stricter when
the number of equations increases.

If instead the exact aligned run identifies a single row near 1e-6, this
protocol is not used; that equation path must be repaired.

## Reference normalization

The accepted H02 B/C screening grid is 20x8x1. The fully compositional two-
component formulation has 8 equations per cell, hence

```
N_ref = 20 * 8 * 8 = 1280 equations
```

The historical strict absolute SNES gate was

```
||R_scaled||_2 <= 1e-6
```

on that grid. Its exactly equivalent per-equation RMS threshold is

```
RMS_ATOL = 1e-6 / sqrt(1280)
         = 2.795084971874737e-8
```

This value is derived only from the already registered 20x8 tolerance and
equation count. It is not fitted to the 60x20 failure.

## New absolute convergence gate

Absolute convergence now requires ALL THREE conditions:

```
||R_scaled||_2 / sqrt(N_global) <= 2.795084971874737e-8
||R_scaled||_infinity          <= 1.0e-6
abs(sum_cells R_mass,i_raw)    <= 8.0e-12 kg/s   for every conserved component i
```

The infinity bound preserves the strongest row-wise implication of the old
global L2 gate: no individual scaled equation may exceed 1e-6.

The component-wise signed mass gate was added after run #86 showed that RMS+Linf
alone completed 60x20 B/C but accumulated Heavy trajectory errors of
approximately 3.13e-6 and 3.28e-6 relative. The 8e-12 kg/s threshold was not
fitted to those failures. It is derived from the already registered 1e-6
trajectory mass budget, the conservative ~0.055 kg conserved-mass scale and
the fixed ~6000 s two-PVI horizon:

```
0.055 kg * 1e-6 / 6000 s ~= 9.2e-12 kg/s
```

and rounded downward to 8e-12 kg/s to cover the 6060-s output overshoot.

The signed sum is formed from the **unscaled production mass-balance rows** after
undoing the solver row scaling and then MPI-summing over owned cells. Internal
face fluxes cancel globally. In fully compositional B/C, H2O and OIL_HEAVY are
checked independently. In A, OIL_HEAVY and the legacy independent-water
conservation row are checked independently.

The existing relative criterion remains `SNES_RTOL=1e-8`. Step-norm
convergence remains disabled with `SNES_STOL=1e-100`.

For fully compositional B/C, the 20x8 criterion is mathematically identical to
the historical absolute L2 criterion. On 60x20 (9600 equations), the equivalent
global L2 threshold is approximately 2.7386e-6, while every row must still be
below 1e-6.

The same RMS and Linf thresholds are used for A so A/B/C share a common
per-equation acceptance standard. This makes the 20x8 A gate slightly stricter
than its historical 1e-6 global-L2 gate; A must therefore be revalidated on
20x8. The already completed 60x20 A run used a stricter global-L2 threshold and
therefore does not need to be rerun if its saved residual/mass audit remains
valid.

## Required validation after activation

1. Re-run 20x8 A/B/C through >2 actual PVI under the mesh-normalized gate.
2. Require existing H2O/Heavy global mass gates <=1e-6.
3. Compare 1 and 2 PVI RF_H, cumulative Heavy, pressure difference and
   B/C viscosity response with the accepted run49/run64 baseline.
4. Run 60x20 B. Only after B completes, run 60x20 C. Do not rerun the already
   certified 60x20 A trajectory.
5. Require the existing <=1e-6 trajectory mass audit before any grid claim.
6. Apply the already preregistered 20x8->60x20 matched-PVI grid gates and record
   PASS/FAIL without changing thresholds.
7. Run the 36x12 aligned-well A/B/C ensemble regardless of the 20x8->60x20
   grid-gate outcome.
8. Emphasize 36x12->60x20 as the final spatial convergence pair. Use the
   versioned certified run64 60x20 A reference rather than rerunning A.

No physical parameter may change during this sequence.

All outputs remain
`CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION`.
