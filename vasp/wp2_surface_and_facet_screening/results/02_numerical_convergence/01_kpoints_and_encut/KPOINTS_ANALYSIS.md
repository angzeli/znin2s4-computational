# WP2 Stage 02 — basal k-point and cutoff analysis

Reviewed on 15 September 2026 from the completed calculation outputs.

The **10×10×1 → 12×12×1 comparison at 500 eV supports retaining 10×10×1 as the provisional mesh for the next cutoff test** on this fixed, relaxed β(001) slab. The slab-energy change per area is 0.098 meV/Å², and both face-specific work-function changes are below 0.025 eV. These satisfy the corresponding working numerical tolerances in the [WP2 roadmap](../../../ROADMAP.md).

This is a limited k-point assessment. There is **no 600 eV calculation in this set**, no validated matched bulk reference for surface excess energies, and no resolved semiconductor gap. The full Stage 02 acceptance decision remains open.

## What was compared

The two source folders are [beta001/500eV_k10](../../../calculation/02_numerical_convergence/01_kpoints_and_encut/beta001_kpoints/500eV_k10) and [beta001/500eV_k12](../../../calculation/02_numerical_convergence/01_kpoints_and_encut/beta001_kpoints/500eV_k12). Energies, forces and completion evidence were read from `OUTCAR`, `OSZICAR` and `vasprun.xml`; eigenvalues and occupations from `EIGENVAL` and XML; vacuum potentials and densities from `LOCPOT` and `CHGCAR`.

Both calculations use the [A2 final geometry](../../../calculation/01_runtime_pilot/01_beta_basal/02_relaxation/CONTCAR): one seven-atom ZnIn₂S₄ septuple layer, surface identity `beta_001_t07`, with its lower and upper faces retained. A common rigid translation of +0.060177578 Å along z centers the slab; the cell, atom order and internal geometry are preserved.

| Common geometry | Value |
| --- | ---: |
| In-plane area, A | 13.02596522 Å² |
| Cell height | 29.41985821 Å |
| Outer-atom slab thickness | 9.45610547 Å |
| Atom-free periodic-image separation | 19.96375274 Å |

The POSCAR geometries are identical, the local POTCAR files are byte-identical, and parsed INCAR settings differ only in `SYSTEM`. Both use VASP 6.6.1, PBE+D3(BJ), the same S/In_d/Zn PAW datasets, ENCUT = 500 eV, and fixed-geometry static calculations. Relevant common controls are `EDIFF=1e-6`, `ISPIN=1`, Gaussian smearing with `SIGMA=0.05 eV`, and a z-directed dipole correction with `DIPOL=(0.5,0.5,0.5)`. Both start from atomic density and output the Hartree-plus-ionic potential (`LVHAR`).

| Sampling and execution | 500eV_k10 | 500eV_k12 |
| --- | ---: | ---: |
| Unshifted Γ-centered mesh | 10×10×1 | 12×12×1 |
| In-plane reciprocal spacing, approximately | 0.18707 Å⁻¹ | 0.15589 Å⁻¹ |
| Actual weighted k points in outputs | 52 | 74 |
| Electronic iterations | 22 | 24 |
| Normal termination / EDIFF reached | Yes / Yes | Yes / Yes |
| VASP elapsed time | 88.755 s | 132.791 s |

Reciprocal spacing is `|b_i|/N_i`, with reciprocal vectors including 2π. Each XML file is complete and contains one static step; neither run reaches `NELM=150`. Both have empty stderr and no identified fatal-error markers. POSCAR, CONTCAR, XML initial/final structures, LOCPOT and CHGCAR agree geometrically within each run. The denser calculation took 1.50 times as long in these particular runs; this is an observed cost comparison, not a general timing benchmark.

## Energies and residual forces

All differences below are **k12 − k10**. E₀ is VASP's reported `energy(sigma→0)` estimate; F is the finite-smearing free energy (`TOTEN`). OUTCAR and XML agree on both energies to the reported precision.

| Quantity | 500eV_k10 | 500eV_k12 | Change |
| --- | ---: | ---: | ---: |
| E₀ (eV per slab cell) | −28.28071051 | −28.27943402 | +0.00127649 |
| F (eV per slab cell) | −28.28099652 | −28.27981974 | +0.00117678 |
| Maximum atomic force norm (eV/Å) | 0.01924870 | 0.02002467 | +0.00077597 |

The energy changes are **0.097996 meV/Å² for E₀** and **0.090341 meV/Å² for F**, about one tenth of the roadmap's 1 meV/Å² working tolerance. The largest change in an individual atomic force vector is 0.00311472 eV/Å, on atom 2 (S, POSCAR order).

The k12 maximum force is only 0.00002467 eV/Å above the Stage 01 relaxation target of 0.02 eV/Å. Because these are static calculations, this is a residual-force diagnostic, not a failed relaxation. The geometry should be checked against the force target again at the eventual selected settings; it should stay frozen during the present numerical comparisons.

The magnitude of the printed entropy contribution per atom is 0.08172 meV at k10 and 0.11021 meV at k12. These are small, but one common Gaussian width does not establish smearing convergence. VASP's extrapolated E₀ needs a width-convergence check for precise zero-smearing claims, and its forces are consistent with F. See the [VASP smearing guidance](https://vasp.at/wiki/Smearing_technique).

For the roadmap's paired excess quantity,

\[
\Gamma_{\mathrm{pair}}=(E_{\mathrm{slab}}-n e_{\mathrm{bulk}})/A,
\]

the reported ΔE/A would equal ΔΓ only if the same bulk-reference contribution were held fixed. A matched bulk-sampling comparison would instead include `−nΔe_bulk/A`, which has not been evaluated here. Therefore the small slab-energy difference is useful screening evidence, **not a validated surface-excess-energy convergence result**. No absolute surface energy or individual-face surface energy is assigned, and the area normalization is A, not 2A.

## Vacuum-referenced electronic quantities

The analysis averages LOCPOT parallel to the surface and treats the lower and upper vacuum regions separately. For each run, CHGCAR integrates to 61.99999960 electrons, consistent with `NELECT=62`; the matched grids are 56×56×432. Density is normalized in electrons/Å³.

The same distant windows are used in both cases: **z = 1.0215–3.9499 Å** below the slab and **25.4700–28.3983 Å** above it. They exclude 6 Å next to each outermost atom and 1 Å next to the periodic boundary, leaving 44 grid points and 2.9284 Å of width per face.

| Window check | k10 lower | k10 upper | k12 lower | k12 upper |
| --- | ---: | ---: | ---: | ---: |
| Potential slope (eV/Å) | −1.124×10⁻⁴ | +1.460×10⁻⁴ | −1.311×10⁻⁴ | +9.705×10⁻⁵ |
| Potential range (eV) | 0.000376 | 0.000480 | 0.000415 | 0.000327 |
| Maximum absolute density (e/Å³) | 8.02×10⁻⁷ | 6.45×10⁻⁷ | 4.67×10⁻⁷ | 2.68×10⁻⁷ |

All four windows satisfy the [existing electrostatics criteria](../../../script/wp2/electrostatics.py): width ≥2 Å, absolute slope ≤0.005 eV/Å, detrended maximum residual ≤0.02 eV, potential range ≤0.03 eV and maximum absolute density ≤10⁻⁵ e/Å³. Using a 5 Å instead of 6 Å setback also passes on both faces and shifts each vacuum mean by less than 0.00012 eV. Windows beginning only 2 Å from the atoms include density tails; at 4 Å the upper window still marginally exceeds the density threshold. Those wider windows are not used as depleted-vacuum plateaus.

The work-function diagnostic is Φ = V_vac − E_F, computed within each run before comparing meshes. This removes an arbitrary additive potential offset. The use of a field-free, charge-depleted region and the Hartree-plus-ionic potential follows the [VASP work-function guidance](https://vasp.at/wiki/Computing_the_Workfunction).

| Referenced quantity (eV) | 500eV_k10 | 500eV_k12 | Change |
| --- | ---: | ---: | ---: |
| Lower-face Φ | 3.995004 | 4.019505 | +0.024501 |
| Upper-face Φ | 6.466849 | 6.484327 | +0.017478 |
| Upper-minus-lower vacuum step | 2.471845 | 2.464822 | −0.007023 |
| Maximum band-31 energy relative to lower vacuum | −3.994872 | −4.000870 | −0.005998 |
| Maximum band-31 energy relative to upper vacuum | −6.466717 | −6.465692 | +0.001025 |

Both Φ changes are below the roadmap's 0.05 eV working tolerance, and the substantial face asymmetry persists. The band-32 minima are only 0.000010 eV above the tabulated band-31 maxima in each case and therefore have essentially identical shifts. These are slab frontier-state diagnostics, not established bulk-like VBM/CBM values. Flat distant vacuum at this cell size does not establish convergence with respect to vacuum length or inter-slab dispersion interactions.

## Frontier bands and the scientific interpretation

With 62 electrons and `ISPIN=1`, the nominal filled-band count is 31. To avoid discarding partially occupied states, the sampled separation is calculated as `min(E_band32) − max(E_band31)` across all sampled k points.

| Frontier diagnostic | 500eV_k10 | 500eV_k12 |
| --- | ---: | ---: |
| Sampled band-31/32 separation | ≈0.000010 eV | ≈0.000010 eV |
| Location of both extrema | Γ | Γ |
| Γ occupations, bands 31 / 32 (full = 1) | 0.498514 / 0.498402 | 0.299067 / 0.298968 |
| Sampled entries with 0.001 < occupation < 0.999 | 2 | 8 |

The approximately 10 µeV separation is not a scientifically resolved gap, especially with 0.05 eV smearing. Both meshes support the same description: **frontier bands touch within numerical resolution and carry partial occupations**. The occupation distribution changes with sampling; the number of partially occupied entries also reflects the different grids and is not a count of physical carriers. Eigenvalues and occupations in EIGENVAL and XML agree within their different printed precision (maximum absolute difference 5×10⁻⁵).

A gap obtained by excluding fractional occupations would be misleading. These results do not establish a semiconductor gap, bulk band alignment, or redox suitability. The single septuple layer has no demonstrated bulk-like interior. Surface-state localization has not been validated here, spin polarization was not tested, and frozen-geometry calculations cannot establish reconstruction stability.

## Decision and next comparison

| Decision | Evidence-based status |
| --- | --- |
| Completion and comparability | Both statics are complete, electronically converged and controlled for geometry and method. |
| k-point sensitivity of the measured slab diagnostics | Within the working energy-per-area and face-specific electronic tolerances for this 10→12 comparison at 500 eV. |
| Provisional mesh | Retain 10×10×1 for the next fixed-geometry comparison. A 16×16×1 test is not compelled by the measured changes, but remains appropriate if tighter near-Fermi or energy claims are needed. |
| Cutoff convergence | Untested. Next compare 500 versus 600 eV at the same 10×10×1 mesh and unchanged geometry, smearing and dipole settings. |
| Full Stage 02 acceptance | Pending matched bulk-reference evidence and the remaining cutoff, vacuum, thickness and final-geometry checks. Spin/localization claims require their own evidence. |

The present recommendation applies to this β(001) model. It does not establish settings for the edge slab or other terminations. No calculations were prepared, launched or released as part of this analysis.
