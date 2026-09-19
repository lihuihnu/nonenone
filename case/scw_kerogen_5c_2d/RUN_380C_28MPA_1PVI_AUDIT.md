# Five-component PR flow screening audit: 380 C / 28 MPa / 60x20 / 1 PVI

Audit date: 2026-09-19.

## Scope

This audit covers GitHub Actions run `35434840156`, the first native 2-D
flow screening of the current experiment-driven five-component topology:

1. H2O
2. OIL_GASOLINE
3. OIL_DIESEL
4. OIL_MIDDLE
5. OIL_HEAVY

The run is explicitly:

`CONDITIONAL_5C_FLOW_SCREENING_NOT_FORMAL_VALIDATION`.

It is not the formal 360/380 C paired laboratory claim and does not promote the
current H2O-lump PR screening BIPs or transport closure to validated status.

## Provenance

- workflow run: `35434840156`
- workflow conclusion: `success`
- artifact: `10583086449`
- artifact SHA-256:
  `7e3471c06456533b0d0dbee76c1b8129124bd7349f656606ee2a3fae805b6f0b`
- head:
  `04d783543a55f5be1c747b65f05a3796d5f82863`

All workflow stages passed:

- case build;
- 20x8 short onset / producer-output-contract smoke;
- 60x20 PR trajectory to 1 actual PVI;
- five-component mass/selectivity audit;
- artifact upload.

## Baseline

- T = 653.15 K (380 C);
- pressure center = 28 MPa;
- PR76;
- 60x20x1 grid over 0.30 x 0.10 x 0.010 m;
- phi = 0.25;
- k = 1e-12 m2;
- effective PV = 7.5e-5 m3;
- pure-H2O injector target = 2.5e-8 m3/s reservoir volume;
- injector max BHP = 30 MPa;
- producer BHP = 28 MPa;
- initial z_H2O = 0.20;
- remaining 0.80 uses the characterized recovered-oil mole ratio.

The resulting initial hydrocarbon mass fractions reproduce the registered
experimental recovered-oil basis exactly to the reported precision:

| component | initial HC mass fraction |
|---|---:|
| Gasoline | 0.008100 |
| Diesel | 0.237300 |
| Middle | 0.341100 |
| Heavy | 0.413500 |

## Numerical health

Final actual PVI:

`1.0000000000279312`.

Solver summary:

- accepted internal steps: 1556;
- rejected/retried steps: 35;
- minimum accepted dt: about 0.0469 s;
- maximum accepted dt: 2.0 s;
- maximum accepted solve SNES iterations: 11.

No well-control switch occurred.

The injector remained on reservoir-rate control throughout the recorded
trajectory. The maximum injector BHP is the initialized 28.05 MPa value and
the converged injector BHP rapidly falls close to 28 MPa. The 30 MPa cap is
never approached.

At 1 PVI:

- injector BHP ~= 28.001302 MPa;
- producer BHP = 28.000000 MPa;
- DeltaP ~= 0.001302 MPa = 1.302 kPa.

The very small pressure drop is consistent with the current high-temperature
LBC mobility package and must not be interpreted as experimentally validated
until the transport-property gate is closed.

## Strict component mass gate

Registered limit: `1e-6` relative to initial+injected component inventory.

| component | max relative error | status |
|---|---:|---|
| H2O | 7.8474e-8 | PASS |
| OIL_GASOLINE | 1.2138e-7 | PASS |
| OIL_DIESEL | 1.2424e-7 | PASS |
| OIL_MIDDLE | 1.2424e-7 | PASS |
| OIL_HEAVY | 1.2424e-7 | PASS |

Therefore:

`FIVE_COMPONENT_STRICT_MASS_GATE = PASS`.

## Matched-PVI recovery/selectivity

### 0.25 PVI

All four hydrocarbon recovery fractions are effectively identical:

- Gasoline: 0.271050;
- Diesel: 0.271050;
- Middle: 0.271050;
- Heavy: 0.271050.

`E_L/H,instant = 1.000000`.

`E_L/H,cumulative = 1.000000`.

There is no measurable early selectivity.

### 0.50 PVI

Recovery remains essentially proportional:

- Gasoline: 0.536495;
- Diesel: 0.536493;
- Middle: 0.536493;
- Heavy: 0.536493.

Gasoline - Heavy recovery difference:

`+0.00000203` = about `+0.000203 pp`.

Cumulative `E_L/H` remains 1.000000 to the reported precision.

### 0.75 PVI

Selectivity begins to emerge:

- Gasoline RF = 0.695864;
- Diesel RF = 0.690239;
- Middle RF = 0.690241;
- Heavy RF = 0.690237.

Gasoline - Heavy:

`+0.005628` = **+0.563 pp**.

Instantaneous:

`E_L/H = 1.003965`.

Cumulative:

`E_L/H = 1.000273`.

### 1.00 PVI

- Gasoline RF = **0.765388**;
- Diesel RF = 0.749113;
- Middle RF = 0.749118;
- Heavy RF = **0.749106**;
- total hydrocarbon RF = 0.749244.

Gasoline - Heavy:

`+0.016283` = **+1.628 pp**.

Diesel - Heavy:

`+7.27e-6` = essentially zero.

Middle - Heavy:

`+1.23e-5` = essentially zero.

The registered Light group is Gasoline + Diesel. At 1 PVI:

- Light-group RF = 0.749650;
- Heavy-group RF = 0.749106;
- instantaneous `E_L/H = 1.008797`;
- cumulative `E_L/H = 1.000727`.

Thus the instantaneous producer stream is only about **0.88% lighter** on the
registered L/H metric, while the cumulative enrichment is only about
**0.073%**.

The five-component baseline therefore does **not** show a strong whole-Light
selective-extraction signal.

## Producer hydrocarbon composition at 1 PVI

At 1 PVI the total producer stream is water-rich:

- Y_H2O = 0.716295.

After renormalizing only the four hydrocarbon components, the instantaneous
hydrocarbon product is:

| component | initial oil mass fraction | producer HC-normalized mass fraction | relative enrichment |
|---|---:|---:|---:|
| Gasoline | 0.008100 | 0.010208 | **1.2603** |
| Diesel | 0.237300 | 0.236800 | 0.9979 |
| Middle | 0.341100 | 0.340410 | 0.9980 |
| Heavy | 0.413500 | 0.412582 | 0.9978 |

So Gasoline itself is about **26% enriched within the instantaneous hydrocarbon
portion** at 1 PVI.

This strong relative enrichment of a very small initial 0.81 wt% cut has only
a small effect on the aggregate Light/Heavy metric because Diesel dominates
the registered Light group.

## Phase-resolved transport mechanism

The phase-resolved producer ledger shows a large ordering in direct
Water-phase carriage.

At 1 PVI, fraction of each cumulative produced hydrocarbon component assigned
to the Water phase:

| component | Water-phase share of cumulative production |
|---|---:|
| Gasoline | **1.989%** |
| Diesel | 9.09e-6 |
| Middle | 1.54e-5 |
| Heavy | 1.61e-8 |

Therefore the observed Gasoline selectivity is not merely a normalization
artifact.

At the 1-PVI producer completion cell the PR water/oil partition coefficients
are approximately:

| component | K_water/oil |
|---|---:|
| H2O | 1.6617 |
| Gasoline | **1.2707e-2** |
| Diesel | 5.72e-6 |
| Middle | 9.69e-6 |
| Heavy | **1.01e-8** |

The full-field median Gasoline `K_water/oil` is about `1.27e-2`, while the
Heavy median is about `1.01e-8`.

This provides a direct conditional mechanism:

`Gasoline has a much larger aqueous/oil partition tendency than the heavier
registered lumps -> a small but finite fraction enters the Water phase ->
Gasoline is produced slightly earlier than Diesel/Middle/Heavy.`

## Phase topology

- initial state: 100% Oil-only cells;
- 0.5 PVI: about 76.58% Oil+Water cells and 23.42% Oil-only cells;
- 1.0 PVI: about 99.17% Oil+Water cells and 0.83% Oil-only cells;
- no active Gas phase is present in the recorded 1-PVI field.

Therefore the Gasoline selectivity is **not a gas-phase liberation effect** in
this screening run. It is primarily associated with water/oil partitioning.

## Interpretation

This first five-component native trajectory establishes that:

1. the experiment-driven five-component production kernel can initialize and
   run through 1 actual PVI in native 2-D flow;
2. first-class producer component output works for all registered oil lumps;
3. five independent component mass ledgers close below the registered limit;
4. component-selective transport is resolvable by the current output contract;
5. the present PR screening baseline predicts a small, physically traceable
   preference for the Gasoline cut after substantial water production begins.

It does **not** establish that real Chang-7 oil will show a +1.63-pp Gasoline
RF advantage or a +26% instantaneous Gasoline enrichment. The H2O-Gasoline,
H2O-Diesel, H2O-Middle and especially H2O-Heavy interaction parameters remain
screening/proxy inputs, and the high-temperature transport closure is not yet
experimentally validated.

## Next scientific step

The next useful numerical experiment is the matched **360 C / 28 MPa**
five-component control under the same geometry, initial composition, PVI rate,
well controls and relative-permeability package.

For that comparison, the PR H2O-lump BIPs must use the already registered
temperature dependence

`kij(T) = kref + b * (1/T - 1/653.15 K)`

rather than freezing the 380 C values.

Then the paired analysis can report the registered formal quantities:

- `Delta RF_i = RF_i,380 - RF_i,360`;
- `E_L/H,380 / E_L/H,360`;
- pressure-drop difference;
- injector-BHP difference;
- saturation/front differences;
- independent component mass closure.

That pair would still be a **conditional mechanism-isolation numerical pair**,
not the formal laboratory validation, until apparatus/specimen and property
calibration gates close.
