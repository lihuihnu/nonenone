# Heringer 2025 water/CO2/BSB-oil three-phase flash

This offline benchmark reproduces the full three-phase PR flash reported at
650 K and 390 bar in Table 6 of:

- J. Heringer et al., *Fluid Phase Equilibria* 594 (2025) 114378,
  <https://doi.org/10.1016/j.fluid.2025.114378>.

The feed and PR component parameters are transcribed from Table B3, all
non-zero classical one-fluid BIPs from Table B4, and the gas/oil/water phase
compositions and mole fractions from Table 6. The extracted values are kept in
`reference_data/` exactly at the printed precision.

The attached paper does not list critical volumes or molecular weights because
they do not enter its PR fugacity calculation. The executable supplies the BSB
values tabulated by Fernandes, Marcondes and Sepehrnoori, *Journal of
Computational Physics* 435 (2021) 110263,
<https://doi.org/10.1016/j.jcp.2021.110263>. Water values are standard pure
component metadata. These fields satisfy `CompositionalMixture` metadata
requirements but do not affect this PR flash result.

The paper cites the original 1976 Peng-Robinson equation and does not state
whether the later high-acentric-factor alpha correction was used. The benchmark
therefore evaluates both PR76 and PR78 without altering the published BIPs. For
each alpha variant it runs the unrestricted production flash and a production
three-phase restricted flash initialized with the rounded Table 6 phase
compositions. The latter distinguishes initialization/basin effects from EOS
parameter differences.

The same state is also calculated with two independent production backends:

- SW1992 uses zero salinity, the original aqueous CO2-H2O correlation and the
  published hydrocarbon correlation for C1 and C2-3. Because that correlation
  was fitted only from methane through n-butane, C4-6 and heavier aqueous BIPs
  remain at the paper's explicit 0.5000 value. A clearly labeled sensitivity
  case extrapolates the correlation to all heavy pseudo-components. Two
  additional attribution cases retain the paper's hydrocarbon BIPs while
  enabling either the SW water alpha plus CO2 correlation, or only the SW water
  alpha, so a loss of the reported three-phase state can be traced to a
  specific SW parameter substitution rather than mistaken for solver failure.
- The principal internally consistent CPA case is the PR-CPA model of
  Sorensen et al., *Journal of Natural Gas Engineering* 3 (2018), 1-38,
  <https://doi.org/10.7569/JNGE.2018.692501>. It uses their jointly regressed
  4C water parameters, simplified radial distribution, temperature-dependent
  water/non-aqueous BIPs, and their one-acceptor-site CO2 solvation model. The
  C2-3 pseudo-component uses the mean of the published C2 and C3 BIPs; C4-6 is
  mapped to nC5; all heavier fractions use the published C7+ value. Their pure
  water calibration reaches 643 K, but the binary mutual-solubility data stop
  below this 650 K case, so this remains a documented extrapolation.
- A second physically motivated heavy-hydrocarbon check uses canonical
  SRK-CPA 4C water and maps the C7+ BSB fractions by molecular weight to the
  water/n-alkane BIPs fitted directly to L-V-W three-phase curves by Jia and
  Okuno, *AIChE Journal* 64 (2018), 3429-3442,
  <https://doi.org/10.1002/aic.16191>. Their highest three-phase measurements
  reach 644 K, but the light BSB fractions and CO2 remain outside that fit.
- The original SRK-CPA case (standard 4C water with all Table B4 BIPs) and the
  PR78-cubic-plus-association case are retained only as transfer/sensitivity
  calculations. PR-fitted BIPs are not asserted to be transferable CPA
  parameters, and the latter water cubic/association parameters were not
  jointly regressed.

`physical_selection.csv` separates target-specific recommendations from
predictive applicability checks. In particular, a converged one-phase result
from a literature CPA parameter set is retained as a model prediction; it is
not replaced by an imposed three-phase solution merely to resemble Table 6.

## Standard-parameter comparison policy

The primary numerical experiment deliberately performs no fit to Table 6 or
to the PR response. `summary.csv` and `phase_comparison.csv` contain exactly
one unrestricted calculation from each EOS family:

- `PR78` uses the published feed, pure-component properties and Table B4 BIPs.
  This is the reproduction check against the paper's PR result.
- `SW_water_alpha_paper_BIPs` changes only the aqueous water alpha treatment
  to the zero-salinity Soreide-Whitson model. It retains the published Table B4
  BIPs and does not extrapolate the SW hydrocarbon correlation to unsupported
  heavy pseudo-components.
- `SRK-CPA_4C-water` uses the standard SRK physical term and full, unscaled 4C
  water association parameters. Table B4 BIPs are transferred unchanged so
  this remains a controlled EOS comparison, not a CPA-specific BIP regression.

The SW and CPA results are therefore independent model predictions and are not
expected to reproduce Table 6. Additional SW correlation substitutions,
literature CPA mappings and seeded solver paths are written only to
`diagnostic_summary.csv`; they do not enter the primary comparison or figure.
`standard_model_basis.csv` records the parameter provenance and confirms that
all three primary cases have `target_fitted=0`.

Run from `tools/`:

```text
make run-heringer2025-bsb
```

Generated comparison files are written to `results/` and intentionally remain
untracked.

The publication-style comparison figure can be regenerated with:

```text
python example/heringer2025_bsb_three_phase/plot_eos_comparison.py \
  --input-dir example/heringer2025_bsb_three_phase/results \
  --output-dir example/heringer2025_bsb_three_phase/results/figures
```
