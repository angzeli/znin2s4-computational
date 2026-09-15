# WP2 Stage 02 — basal cutoff convergence

Reviewed on 15 September 2026 from the completed [500 eV](../../../calculation/02_numerical_convergence/01_kpoints_and_encut/beta001_encut/500eV_k10) and [600 eV](../../../calculation/02_numerical_convergence/01_kpoints_and_encut/beta001_encut/600eV_k10) calculations.

**Retain 500 eV provisionally for this fixed β(001) model.** Increasing the cutoff to 600 eV changes the slab energy by only 0.049 meV/Å² and both face-specific work-function diagnostics by less than 0.002 eV. Residual forces change little, and the frontier bands remain near-touching. These results support the cutoff choice for the measured slab diagnostics; matched bulk references and the remaining Stage 02 checks are still required.

## Controlled comparison

Both calculations use the same frozen, centered [A2 geometry](../../../calculation/01_runtime_pilot/01_beta_basal/02_relaxation/CONTCAR): one seven-atom ZnIn₂S₄ septuple layer, surface `beta_001_t07`, with lower and upper face identities preserved. The area is 13.02596522 Å², slab thickness 9.45610547 Å and atom-free periodic-image separation 19.96375274 Å. This is the same geometry used in the [k-point comparison](KPOINTS_ANALYSIS.md).

POSCAR, KPOINTS and local POTCAR bytes are identical between the two runs. Parsed INCAR settings differ only in `ENCUT` and the descriptive `SYSTEM` label. Both use VASP 6.6.1, PBE+D3(BJ), a Γ-centered 10×10×1 mesh, `EDIFF=10⁻⁶ eV`, `ISPIN=1`, Gaussian smearing with `SIGMA=0.05 eV`, and z-directed dipole correction centered at `(0.5,0.5,0.5)`. They are fresh fixed-geometry statics (`NSW=0`, `ISTART=0`, `ICHARG=2`).

Both terminated normally, reached EDIFF, and contain one complete static step. XML and OUTCAR agree on energies, forces and executed inputs within printed precision; initial/final structures and density/potential geometries agree. Neither has timeout or identified fatal-error evidence, and both stderr files are empty.

| Execution | 500 eV | 600 eV |
| --- | ---: | ---: |
| Electronic iterations | 22 | 24 |
| Weighted k points | 52 | 52 |
| VASP elapsed time (s) | 107.704 | 138.972 |
| Density/potential grid | 56×56×432 | 64×64×480 |
| Integrated charge (electrons) | 61.99999960 | 62.00000067 |

The denser real-space grid follows the higher cutoff at unchanged `PREC=Accurate`. Density and potential grids match within each run; the analysis uses physical vacuum intervals on each native grid. No interpolation is applied.

## Energetic and electronic sensitivity

All changes below are **600 − 500 eV**. E₀ denotes VASP's `energy(sigma→0)` estimate; F is the finite-smearing free energy (`TOTEN`).

| Quantity | 500 eV | 600 eV | Change |
| --- | ---: | ---: | ---: |
| E₀ (eV/slab cell) | −28.28071051 | −28.28134809 | −0.00063758 |
| F (eV/slab cell) | −28.28099652 | −28.28163407 | −0.00063755 |
| Maximum force (eV/Å) | 0.01924870 | 0.01891137 | −0.00033733 |
| Lower-face Φ (eV) | 3.995004 | 3.993362 | −0.001641 |
| Upper-face Φ (eV) | 6.466849 | 6.465011 | −0.001838 |
| Upper-minus-lower vacuum step (eV) | 2.471845 | 2.471648 | −0.000197 |
| Band-31 maximum relative to lower vacuum (eV) | −3.994872 | −3.993232 | +0.001640 |
| Band-31 maximum relative to upper vacuum (eV) | −6.466717 | −6.464880 | +0.001837 |

The changes in E₀/A and F/A are −0.048947 and −0.048945 meV/Å². The largest change in an individual atomic force vector is 0.00042909 eV/Å, on atom 7 (POSCAR order). Both residual maximum forces are below the 0.02 eV/Å pilot target; these static checks do not constitute new relaxations. The 500 eV repeat reproduces the earlier k10 E₀ to printed precision. Smearing is held fixed, so no zero-smearing convergence claim follows.

Φ is calculated as `V_vac − E_F` within each run before comparing cutoffs. The Hartree-plus-ionic potential (`LVHAR`) is averaged parallel to the slab, with each face treated separately, following the [VASP work-function guidance](https://vasp.at/wiki/Computing_the_Workfunction).

The comparison uses the same 6 Å setback from each outermost atom and 1 Å exclusion from the periodic boundary. Actual sampled intervals are 1.0215–3.9499 / 25.4700–28.3983 Å at 500 eV and 1.0420–3.9226 / 25.4972–28.3779 Å at 600 eV (lower/upper). All four windows pass: absolute slopes ≤0.000146 eV/Å, residuals ≤0.000104 eV, potential ranges ≤0.000481 eV and maximum absolute densities ≤8.02×10⁻⁷ e/Å³. Window widths exceed 2.88 Å.

All 16 windows at 2/4/5/6 Å setbacks were checked. Both faces pass at 5 and 6 Å; changing between these setbacks shifts each vacuum mean by less than 0.00012 eV. At 2 Å, both faces fail residual, range and density criteria; at 4 Å, the upper face still fails the density criterion. No jumps were detected inside the tested intervals. These screening results do not establish vacuum-size convergence.

With 62 electrons, the frontier diagnostic is `min_k E[32,k] − max_k E[31,k]`, retaining partially occupied states. Both cutoffs give approximately +0.000010 eV, with both extrema at Γ. Γ occupations are 0.498514/0.498402 at 500 eV and 0.498527/0.498411 at 600 eV (full occupation = 1). Complete EIGENVAL/XML arrays agree within printed precision, and weighted occupations reproduce NELECT. The separation remains below the 0.0001 eV reporting threshold: **no resolved semiconductor gap is demonstrated**. Band-32 minima lie only 0.000010 eV above the tabulated band-31 maxima; their vacuum-referenced shifts are the same to the shown precision.

## Decision and limits

The measured slab-energy and vacuum-referenced changes are comfortably within the [roadmap's working tolerances](../../../ROADMAP.md) of 1 meV/Å² and 0.05 eV. There is no evidence here requiring 600 eV for these diagnostics on this geometry. Cutoff adequacy remains specific to the quantity being tested, as emphasized in the [VASP ENCUT guidance](https://vasp.at/wiki/ENCUT).

The roadmap's energetic criterion concerns `Γ_pair=(E_slab−n e_bulk)/A`. This comparison establishes only ΔE_slab/A; it lacks the matched `−nΔe_bulk/A` term. It therefore does **not** certify surface-excess-energy convergence or assign individual-face surface energies. Keep 500 eV and 10×10×1 as provisional settings while assessing vacuum, thickness and final-geometry consistency. The single layer has no demonstrated bulk-like interior; spin state, surface-state localization and reconstruction stability were not tested by these frozen, nonmagnetic calculations.
