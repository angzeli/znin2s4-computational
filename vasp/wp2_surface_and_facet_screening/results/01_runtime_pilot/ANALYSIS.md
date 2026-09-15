# Stage 01 — Runtime pilot results

All five pilots terminated normally and converged electronically. **A2 is a force-converged basal seed; B1 requires further relaxation and structural review.** Results refer to the calculations reviewed on 15 September 2026.

| Pilot | Role | Ionic steps / SCF | Final maximum force (eV/Å) | Frontier separation (eV) |
| --- | --- | ---: | ---: | ---: |
| A0 | Uncorrected basal static | 1 / 29 | 0.5113791 | −0.116564 |
| A1 | Dipole-corrected basal static | 1 / 20 | 0.5290504 | −0.127949 |
| B0 | Unrelaxed edge static | 1 / 41 | 2.0424329 | −0.073611 |
| A2 | Relaxed basal seed | 18 / 139 | 0.0198316 | +0.000010 |
| B1 | Unfinished edge relaxation | 100 / 1772 | 0.0374145 | +0.742242 |

## Interpretation

A2 meets the 0.02 eV/Å force target and retains its septuple sheet and cation coordination. Its maximum displacement is 0.190604 Å; final thickness and atom-free separation are 9.456105 and 19.963753 Å. This single-layer seed has no demonstrated bulk-like interior.

B1 reached NSW=100 without ionic convergence or timeout evidence. Its lowest maximum force was 0.0358072 eV/Å at step 99; the endpoint remains step 100. The maximum displacement is 1.679712 Å, and the S–S contact between sites 6 and 8 shortens from 3.817295 to 2.094616 Å. Five cation–S connections are lost and four gained under the original cutoffs. All atoms remain connected, but periodic connectivity depends on the In–S cutoff. Final thickness and atom-free separation are 12.549103 and 19.766104 Å. The full trajectory retains the step-50 force excursion.

Frontier separation is `min_k E[N+1,k] − max_k E[N,k]`, with N=NELECT/2 (31 basal; 124 edge). All sampled states are retained, including partial occupations. EIGENVAL/XML agree within printed precision and their weighted occupations reproduce NELECT using full occupation 1. Negative values indicate sampled overlap; ±0.0001 eV is a reporting threshold for unresolved separation. At A2 Γ, energies are −1.022448/−1.022438 eV and occupations 0.498501/0.498390. B1's minimum direct separation is 0.774108 eV at an unconverged geometry. The former occupation-filtered `sampled_gap_eV` values are superseded. None establishes a production surface band gap or experimental metallicity.

A0 retains a vacuum slope of approximately −0.164 eV/Å. A1/B0/A2/B1 pass distant-window screening, with near-surface warnings retained. Both faces are checked at 2/4/5/6 Å setbacks with 1 Å boundary exclusion: width ≥2 Å, absolute slope ≤0.005 eV/Å, residual ≤0.02 eV, range ≤0.03 eV and maximum absolute density ≤10⁻⁵ e/Å³. Passing windows do not establish vacuum/thickness convergence or equivalent faces.

## Data and figures

- [Pilot summary](runtime_pilot_summary.csv): corrected diagnostics and status, plus displacement, free-energy, entropy and dipole metrics from the independent audit.
- [Vacuum checks](vacuum_window_checks.csv): all 40 windows, including failures.
- [Planar profiles](planar_profiles.csv): all five pilots, retaining signed density values.
- [Relaxation histories](relaxation_history.csv): all 18 A2 and 100 B1 steps. Combined tables identify each run in the `pilot` column.
- Figures: [basal statics](figures/basal_vacuum_review.pdf), [edge static](figures/edge_vacuum_review.pdf), [A1/A2 vacuum](figures/A1_A2_vacuum.pdf), [A2 relaxation](figures/A2_relaxation.pdf), [B1 relaxation](figures/B1_relaxation.pdf).

[Calculation inputs and outputs](../../calculation/01_runtime_pilot/) retain the executed settings and structures. [Regression evidence](../../script/wp2/fixtures/STAGE01_PILOT_EVIDENCE.md) records independent validation. Numerical, vacuum and thickness convergence remain necessary before final surface comparisons.
