# Verified active-root / candidate-Gibbs separation

## Revision and evidence

Source fix: e03c80242b362a994038ceb7df00afd3bf01567f.
Remote verification run: 35325822571, job 105538507913, SUCCESS.
Artifact: 10538574579 (cpa-root-policy-verified).
Artifact ZIP SHA256: cf3fac2dffadcc6df8de788f87753c7d8aae43629899bf85fc1dce28c305ebe1.
The artifact contains the applied patch, source parent, resulting verified commit, branch CSV certificates, and failure-mode regression log.
The temporary write-enabled patch-application workflow was removed after its successful non-force push; normal read-only regression workflows remain.

## Implemented contract

Active Oil and Water phases use liquid-type density roots, and active Gas uses the vapor-type root throughout residual evaluation, Newton iterations, branch continuation, canonicalization, and candidate total-Gibbs comparison. The cpaSelectGibbsMinimumRoot option now affects missing-phase candidate evaluation only. The active reference chemical potentials never silently switch to a different root.

A converged restricted OW branch that fails missing-phase stability is not accepted as equilibrium. For the opt-in CPA recovery path it supplies an initial state for active-set expansion; the expanded state must converge, close equilibrium/material equations, and pass stability before its total Gibbs energy is compared. No unstable restricted state is relabeled as stable.

All Heavy and Athabasca parameter values, EOS formulas, and numerical tolerances are unchanged by this five-file source fix. The option default remains false; the new expansion recovery is CPA/opt-in only.

## Actual regression results

The original five fixtures now pass 5/5 both locally and on GitHub Actions (previously 3/5):

| T (K) | P (MPa) | water mass fraction | cold/continued maximum composition gap |
|---:|---:|---:|---:|
| 548.2 | 6.91 | 0.559 | 1.63858e-12 |
| 573.1 | 9.52 | 0.559 | 4.77396e-15 |
| 573.1 | 9.52 | 0.441 | 3.66374e-15 |
| 583.0 | 11.00 | 0.441 | 1.88675e-10 |
| 603.5 | 15.32 | 0.559 | 1.93401e-13 |

Across the 15 recorded states: maximum material closure 2.23155e-14; maximum fixed-root log-fugacity spread 9.50593e-10; stored/recomputed active-root Z mismatch zero. The regression was strengthened: even metastable OW branches must remain root-consistent and seed-independent, and the accepted global candidate must not have higher total Gibbs energy than the restricted branches.

At 603.5 K / 15.32 MPa / water mass fraction 0.559, the unrestricted solution has three phases and passes its numerical certificate. Both restricted OW paths remain unstable. With the same dimensionless Gibbs reference, the three-phase result is -0.67209962007606705 versus -0.67209722991578091 for cold OW. Thus phase release is retained, not suppressed to make a two-liquid test pass.

The existing thermodynamic_failure_modes_test.cpp suite also reports ALL PASS, including its deterministic PR/SW/CPA robustness cases.

## Limits

This certifies the fixed-root bug repair on the registered regression fixtures, not a mathematical proof of global minimization for all feeds. Full Figure-7 pointwise boundary reproduction, nontrivial incipient-vapor identification, and complete boundary path-independence require their separate verification. Aggregate AAD remains observational; it is never a tuning target. OIL_HEAVY physical validation and formal flow promotion remain BLOCKED.
