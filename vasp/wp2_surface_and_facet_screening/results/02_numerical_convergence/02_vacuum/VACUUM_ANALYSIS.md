# WP2 Stage 02 — basal vacuum convergence

Reviewed on 15 September 2026 from the completed [19.964 Å baseline](../../../calculation/02_numerical_convergence/01_kpoints_and_encut/beta001_encut/500eV_k10), [25 Å](../../../calculation/02_numerical_convergence/02_vacuum/beta001_25A) and [30 Å](../../../calculation/02_numerical_convergence/02_vacuum/beta001_30A) calculations.

**Retain approximately 20 Å of actual atom-free separation for this single-layer β(001) model.** The largest energy change from the baseline is 0.110 meV/Å², and the largest face-specific work-function change is 0.00268 eV. Successive refinements pass the [roadmap's working tolerances](../../../ROADMAP.md) of 1 meV/Å² and 0.05 eV. This establishes vacuum insensitivity at those tolerances for the measured quantities on this fixed model; thickness and final-geometry validation remain open.

## Comparison and completion

Vacuum means the outer-atom slab-to-image separation, not the total cell height. The baseline is **19.96375274 Å**, rather than exactly 20 Å. All cases retain the frozen A2 septuple sheet (`beta_001_t07`), its seven atoms and species order, internal Cartesian geometry, lower/upper face identities, 9.45610547 Å thickness and 13.02596522 Å² area. Only the normal cell vector and a common centering translation change.

The KPOINTS and local POTCAR files are byte-identical. INCAR settings differ only in the descriptive `SYSTEM` label: VASP 6.6.1, PBE+D3(BJ), 500 eV, Γ-centered 10×10×1 sampling, `EDIFF=10⁻⁶ eV`, `ISPIN=1`, Gaussian `SIGMA=0.05 eV`, and z-directed dipole correction with `DIPOL=(0.5,0.5,0.5)`. These are fresh static calculations, using the provisional [k-point](../01_kpoints_and_encut/KPOINTS_ANALYSIS.md) and [cutoff](../01_kpoints_and_encut/ENCUT_ANALYSIS.md) selections.

All three terminated normally, reached EDIFF and contain one complete static step. XML/OUTCAR energies, forces and executed settings agree within printed precision. POSCAR, CONTCAR, XML and density/potential structures agree within each run. Stderr is empty, with no identified fatal-error or timeout evidence.

| Quantity | 19.964 Å | 25 Å | 30 Å |
| --- | ---: | ---: | ---: |
| Cell height (Å) | 29.41985821 | 34.45610547 | 39.45610547 |
| Electronic iterations | 22 | 27 | 25 |
| VASP elapsed time (s) | 107.704 | 163.283 | 177.827 |
| Density/potential grid | 56×56×432 | 56×56×504 | 56×56×576 |
| E₀ (eV/slab cell) | −28.28071051 | −28.27980524 | −28.27927799 |
| F (eV/slab cell) | −28.28099652 | −28.28009128 | −28.27956401 |
| Maximum force (eV/Å) | 0.01924870 | 0.01971234 | 0.01967039 |
| Lower-face Φ (eV) | 3.995004 | 3.992328 | 3.992668 |
| Upper-face Φ (eV) | 6.466849 | 6.464842 | 6.465043 |
| Upper-minus-lower vacuum step (eV) | 2.471845 | 2.472515 | 2.472375 |

E₀ is VASP's `energy(sigma→0)` estimate; F is the finite-smearing free energy. Smearing is held fixed, so this does not establish zero-smearing convergence. Timing differences are observations from these runs, not a general benchmark.

## Changes and energetic meaning

Changes are larger-vacuum minus smaller-vacuum results, normalized by the single in-plane area A.

| Comparison | ΔE₀ (eV) | ΔE₀/A (meV/Å²) | ΔΦ lower (eV) | ΔΦ upper (eV) | Largest atomic force-vector change (eV/Å) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 19.964 → 25 Å | +0.00090527 | +0.069497 | −0.002676 | −0.002007 | 0.00143335 |
| 25 → 30 Å | +0.00052725 | +0.040477 | +0.000340 | +0.000201 | 0.00105281 |
| 19.964 → 30 Å | +0.00143252 | +0.109974 | −0.002336 | −0.001806 | 0.00248616 |

Both successive energy changes and the cumulative change are below 1 meV/Å². The corresponding F changes differ from ΔE₀ by at most 0.00000003 eV. All residual maximum forces remain below the 0.02 eV/Å pilot target; a fixed-geometry force check does not constitute a new relaxation.

For `Γ_pair=(E_slab−n e_bulk)/A`, composition, area, functional, PAW datasets, cutoff and in-plane sampling are unchanged. With **one common bulk reference**, the bulk term cancels, so `ΔΓ_pair=ΔE_slab/A` for this vacuum-only comparison. No absolute Γ value or individual-face surface energy is assigned. This cancellation does not resolve the matched-reference questions in comparisons that vary cutoff or sampling.

Energy was checked independently of vacuum flatness. The printed D3 contribution changes from −1.74582 to −1.74525 to −1.74501 eV; the total energetic assessment includes this residual vacuum dependence. A flat potential alone would not establish that dependence is negligible.

## Vacuum reference and frontier states

Φ is evaluated as `V_vac−E_F` within each run, using the plane-averaged Hartree-plus-ionic potential (`LVHAR`) and treating the two faces separately. This follows the [VASP work-function guidance](https://vasp.at/wiki/Computing_the_Workfunction). Raw Fermi levels shift substantially with cell size (−1.019431, −1.633275 and −2.089076 eV); those unaligned values are not used to infer band shifts.

The reported references use a 6 Å setback from each outermost atom and 1 Å exclusion from each periodic boundary. Actual window widths are 2.928, 5.469 and 7.946 Å per face. Native grids are used without interpolation; matching CHGCAR grids integrate to 62 electrons within 1.4×10⁻⁶ electrons. Signed density oscillations are retained.

All six distant windows pass the established screening criteria. Their absolute slopes are ≤0.000146 eV/Å, detrended residuals ≤0.000240 eV, potential ranges ≤0.000481 eV and maximum absolute densities ≤8.02×10⁻⁷ e/Å³. The respective screening limits are 0.005 eV/Å, 0.02 eV, 0.03 eV and 10⁻⁵ e/Å³, with width ≥2 Å.

All 24 windows at 2/4/5/6 Å setbacks were checked. The pattern is unchanged across the three cells: both faces have density warnings at 2 Å, the upper face retains a density warning at 4 Å, and both faces pass at 5 and 6 Å. No discontinuities occur inside the tested intervals. Changing from 6 to 5 Å setback shifts each vacuum mean by less than 0.00012 eV. The substantial face asymmetry persists; the two vacuum levels are not averaged into a single reference.

At all three separations, `min_k E[32,k]−max_k E[31,k]` is approximately +0.000010 eV, with both extrema at Γ and partially occupied frontier states. Complete EIGENVAL/XML arrays agree within printed precision and weighted occupations reproduce NELECT. This remains below the 0.0001 eV reporting threshold: no resolved semiconductor gap is demonstrated. The largest change in a vacuum-referenced band-31 maximum is 0.002677 eV; band-32 minima differ by only 0.000010 eV within each run and show the same sensitivity.

The tested vacuum range is sufficient for the current numerical tolerances on this frozen basal seed. Carry approximately 20 Å forward provisionally, then verify vacuum consistency for the eventual thicker or reconstructed model. Slab thickness, bulk-like interior, spin state, localization and relaxed structural stability remain separate questions.
