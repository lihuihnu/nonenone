# H02 mechanism decomposition: B-A versus C-B

Audit date: 2026-09-19.

## Scope

This note returns from numerical grid certification to the scientific H02
question:

- what does `B-A` represent quantitatively?
- what does `C-B` represent quantitatively?
- where in the 2-D domain do these increments arise?

The primary analysis uses the common certified 60x20 grid so that A/B/C
differences are not mixed with spatial-grid differences.

Sources:

- A: certified run #64 60x20 reference;
- B/C: run #145 60x20 composite-callback revalidation.

All results remain:

`CONDITIONAL_MECHANISM_EXPERIMENT_NOT_REAL_HEAVY_VALIDATION`.

No statement here validates the provisional Heavy identity, H2O-Heavy
`kij=0.30`, or the viscosity laws against real Heavy experiments.

## Experimental counterfactual definitions

- **A**: independent-water / no-transfer counterfactual.
- **B**: fully compositional PR equilibrium enabled; the H02 B viscosity closure
  remains at the registered high-viscosity reference branch for the realized
  oil compositions.
- **C**: same PR equilibrium as B; only the composition-viscosity feedback is
  changed to the current log-mixing law.

Therefore the ordered decomposition is:

- `B-A`: interphase-partitioning / phase-repartitioning increment;
- `C-B`: composition-viscosity / mobility-feedback increment.

This is an ordered counterfactual attribution, not a claim that arbitrary
nonlinear mechanisms possess a unique additive decomposition.

## 1. Quantitative recovery contribution

### Heavy recovery fraction

| PVI | A RF_H | B RF_H | C RF_H | B-A [pp] | C-B [pp] | C-A [pp] |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.309155 | 0.377088 | 0.380319 | +6.7933 | +0.3231 | +7.1164 |
| 1.0 | 0.392886 | 0.463708 | 0.471546 | +7.0823 | +0.7838 | +7.8660 |
| 1.5 | 0.442832 | 0.511501 | 0.521548 | +6.8670 | +1.0046 | +7.8716 |
| 2.0 | 0.479186 | 0.545389 | 0.556994 | +6.6203 | +1.1605 | +7.7809 |

At 1 PVI, the ordered share of the total C-A RF increment is:

- B-A: **90.0%**
- C-B: **10.0%**

At 2 PVI:

- B-A: **85.1%**
- C-B: **14.9%**

Thus the phase-partitioning increment is the dominant H02 recovery mechanism,
while the viscosity feedback is a smaller but progressively accumulating
secondary mechanism.

The B-A RF increment reaches a broad maximum of about **+7.10 pp** near
0.86 PVI and then decreases modestly.

The C-B increment is slightly negative during the earliest transient
(minimum about -0.18 pp near 0.22 PVI), crosses positive at about **0.29 PVI**,
and then grows monotonically to +1.16 pp at 2 PVI. The viscosity benefit
therefore requires a developed mixed/swept region before it appears as a
cumulative Heavy-recovery gain.

## 2. B-A is not direct aqueous Heavy extraction

The phrase "water carries Heavy" would be misleading for the current H02
solution.

At 1 PVI in B:

- total cumulative Heavy production:
  `2.6716e-2 kg`;
- cumulative Heavy assigned to the Water phase:
  `2.7853e-10 kg`;
- Water-phase share of produced Heavy:
  `1.04e-8` as a fraction, or about `1.04e-6 %`.

At 2 PVI the corresponding fraction remains only about `2.31e-8`.

The local PR water/oil Heavy partition coefficient is also tiny. At 1 PVI,
the maximum reconstructed fraction of local Heavy moles in the Water phase is
only O(`1e-7`).

Therefore:

`B-A != direct Heavy dissolved in the aqueous phase and advected out`.

The correct interpretation is **interphase partitioning and phase-volume /
mobility redistribution**, dominated by H2O entering the oil phase while Heavy
remains overwhelmingly in the oil phase.

## 3. What B changes physically

### H2O partitions into the oil phase

At 1 PVI, among cells containing non-negligible water:

- median fraction of local H2O moles residing in the oil phase: about **10.9%**;
- mean: about **16.0%**;
- 90th percentile: about **25.5%**.

At the producer by 1 PVI, about `4.60e-4 kg` of cumulative H2O has been
produced through the Oil-phase ledger, while the Water-phase H2O ledger
contains about `1.544e-2 kg`.

So the important cross-phase transfer is H2O -> oil, not Heavy -> water.

### Saturation redistribution

At 1 PVI, relative to A, B has area-average:

- water saturation change:
  `mean(Sw_B-Sw_A) = -0.0179`;
- oil saturation change:
  `mean(So_B-So_A) = +0.0179`.

Because H02 uses `kr = S^2`, this raises the area-mean oil relative mobility
(`kr_o/mu_o`) by about **7.5%** at 1 PVI.

The B oil viscosity remains exactly the reference `0.002 Pa s` across the
realized field: the oil-phase H2O mass fraction is far below the B closure's
25% activation threshold. Hence this B-A oil-mobility change is a clean
saturation/phase-partitioning effect, not a hidden viscosity reduction.

### Why pressure drop increases in B

At 1 PVI:

- A DeltaP: about `0.003476 MPa`;
- B DeltaP: about `0.003659 MPa`;
- B is about **5.26% higher**.

This is consistent with the phase redistribution.

Although B increases oil-phase mobility, it lowers water-phase saturation.
Because the Water phase has viscosity `5e-5 Pa s`, forty times lower than the
reference oil viscosity, reducing Water-phase relative mobility can outweigh
the increase in oil mobility in the total conductance.

At 1 PVI, an area-average `kr/mu` proxy gives approximately:

- oil mobility: **+7.5%** in B versus A;
- water mobility: **-3.9%**;
- total oil+water mobility proxy: **-3.3%**.

So B produces more Heavy by keeping more of the flowing inventory in the
Heavy-bearing oil phase, not by globally lowering flow resistance.

## 4. Where B-A occurs

At 1 PVI the centroid of the absolute B-A water-saturation difference is
approximately:

`(x,y) = (0.140, 0.078) m`

in the 0.30 m x 0.10 m slab.

Only about **6%** of the absolute B-A saturation difference lies within 3 cm of
the producer. Therefore B-A is not a completion-cell artifact.

The spatial distribution is strongly biased toward the upper swept/front
region:

- lower third of the slab height: about 11% of |Delta Sw|;
- middle third: about 9%;
- upper third: about **80%**.

The positive B-A oil-mobility increment is even more concentrated in the upper
third at 1 PVI: about **88%** of the positive increment.

Thus the phase-partitioning mechanism acts primarily along the upper
water/oil displacement and phase-onset band that connects the swept interior
toward the producer, rather than inside the producer cell itself.

## 5. What C adds: viscosity/mobility feedback

B and C use the same PR equilibrium. Their difference is the viscosity law.

For the realized oil composition, the oil-phase water mass fraction is about
3.9% over most oil-containing cells.

- B: this is below its 25% activation threshold, so
  `mu_o ~= 0.002 Pa s`.
- C: log mixing responds immediately, giving
  `mu_o ~= 0.00173245 Pa s`.

The characteristic C-B oil-viscosity reduction is therefore:

**about 13.38%**.

At 1 PVI, C also has slightly lower oil saturation than B:

- `mean(So_C-So_B) ~= -0.0102`;
- the associated area-mean `kr_o` change is about -2.4%.

Nevertheless the viscosity decrease dominates. The reconstructed oil mobility
`kr_o/mu_o` is at 1 PVI:

- median C/B ratio: about **1.098**;
- area-mean C/B ratio: about **1.102**;
- 90th percentile: about **1.155**.

So the net oil-mobility gain is about 10% despite the slight reduction in oil
saturation.

The area-average total oil+water mobility proxy is about **6.6% higher** in C
than B at 1 PVI.

This appears directly in pressure response:

- B DeltaP: `0.003659 MPa`;
- C DeltaP: `0.003419 MPa`;
- C-B change: **-6.55%**.

At 2 PVI the C-B pressure-drop reduction remains about **-5.99%**.

This is the cleanest numerical signature of the intended "viscosity reduction"
mechanism in H02.

## 6. Where C-B occurs

At 1 PVI the centroid of the absolute C-B saturation redistribution is:

`(x,y) ~= (0.146, 0.045) m`

near the middle of the slab rather than at either well.

Only about **3.4%** of the absolute C-B saturation difference lies within
3 cm of the producer.

Unlike B-A, the saturation redistribution is broadly distributed across the
interior displacement corridor:

- lower height third: about 41%;
- middle third: about 32%;
- upper third: about 26%.

The positive oil-mobility increment itself is more strongly weighted toward
the swept upper/right pathway to the producer; at 1 PVI about 60% of the
positive oil-mobility increment lies in the upper height third and about 41%
in the downstream x-third.

Therefore C-B should be interpreted as a distributed mobility-feedback
mechanism acting throughout the mixed oil/water swept region, with the
strongest transport benefit along the developed high-flow path toward the
producer.

## 7. Mechanistic summary

The current H02 result supports the following conditional causal picture.

### B-A: dominant, thermodynamic/phase-redistribution mechanism

`H2O-oil equilibrium partitioning`
-> H2O enters oil phase
-> Water-phase saturation decreases / oil-phase saturation increases
-> Heavy remains almost entirely in oil
-> Heavy-bearing oil relative mobility increases
-> more Heavy is displaced and recovered

This contributes about **7.08 RF percentage points at 1 PVI** and **6.62 pp at
2 PVI**.

Calling this "aqueous Heavy solubilization" or "Heavy carried by the water
phase" would be incorrect for the current simulation.

### C-B: secondary, transport/viscosity mechanism

`same phase equilibrium`
-> composition-sensitive oil viscosity decreases by about 13.4%
-> net oil mobility rises about 10%
-> pressure drop falls about 6%
-> front/flow-path redistribution increases Heavy production

This contributes about **0.78 RF percentage point at 1 PVI** and **1.16 pp at
2 PVI**.

The relative importance of this viscosity mechanism grows with PVI.

## Interpretation boundary

These are conditional mechanism-isolation results under the frozen H02
surrogate Heavy, PR `kij`, relative permeability, and viscosity closures.

The analysis establishes what the current numerical model is doing. It does
not establish that a real Chang-7 Heavy fraction will partition, swell, or
thin by the same quantitative amounts at 380 C.
