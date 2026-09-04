# WP1 Electronic-Structure KPOINTS Selection Justification

## Purpose

Stage 04 uses two distinct reciprocal-space sampling strategies because density-of-states calculations and conventional band-structure calculations answer different questions:

1. **DOS/PDOS calculations** use dense, uniform, unshifted Gamma-centred meshes to sample the full Brillouin zone and improve spectral smoothness.
2. **Band-structure calculations** use explicit high-symmetry Line-mode paths to resolve eigenvalue dispersion along crystallographically defined reciprocal-space directions.

These strategies are complementary and should not be interchanged. A uniform mesh is appropriate for Brillouin-zone integration but does not define a conventional band path; a Line-mode path resolves selected symmetry directions but is not a replacement for uniform DOS integration.

## Basis from Stage 02 K-Point Convergence

The Stage 02 convergence study established an unshifted Gamma-centred production rule of

`Delta k_max = max_i(|b_i| / N_i) <= approximately 0.156 Angstrom^-1`,

where `b_i` are reciprocal-lattice vectors including the `2 pi` convention and `N_i` are the mesh divisions. The directly validated production meshes were `4 x 4 x 4` for spinel and `12 x 12 x 4` for beta. The spinel energy changed by approximately `0.001481 meV/f.u.` from `4 x 4 x 4` to `5 x 5 x 5`; beta changed by approximately `0.102590 meV/f.u.` from `12 x 12 x 4` to `15 x 15 x 5`. Both are far below the `1 meV/f.u.` target.

The cubic spinel and layered beta series therefore validate both relevant reciprocal-space regimes. The other layered phases have similar in-plane lattice constants, and their different c-axis repeats were accounted for through reciprocal-space intervals rather than by blindly copying an integer mesh. This provided a substantial convergence margin and made additional phase-specific k-point tests or further Stage 03 energy-density refinement unnecessary. The completed Stage 03 calculations used `4 x 4 x 4` for spinel and `12 x 12 x 4` for every layered phase; the latter conservatively oversamples the c direction of alpha1, IIa-prime, and IIb relative to the minimum Stage 02 density recommendation.

## DOS/PDOS KPOINTS

The actual Stage 04 DOS/PDOS files use uniform, unshifted Gamma-centred meshes. The density increase below is the ratio of full-grid point counts before symmetry reduction; irreducible k-point counts may scale differently with space-group symmetry. Directional intervals were calculated from each frozen Stage 04 lattice as `|b_i| / N_i`.

| Phase | Static-SCF Mesh | DOS/PDOS Mesh | Approx. Density Increase | Justification |
| --- | --- | --- | --- | --- |
| Spinel | `4 x 4 x 4` | `6 x 6 x 6` | `3.375x`; `Delta k_max = 0.098629 Angstrom^-1` | Increases sampling by 1.5 in every reciprocal direction for smoother three-dimensional DOS/PDOS integration. |
| Alpha1 | `12 x 12 x 4` | `18 x 18 x 6` | `3.375x`; `Delta k_max = 0.103792 Angstrom^-1` | Refines both the layered in-plane sampling and the already conservative c-axis sampling. |
| Beta | `12 x 12 x 4` | `18 x 18 x 6` | `3.375x`; `Delta k_max = 0.103929 Angstrom^-1` | Uses the densest mesh directly tested in the Stage 02 beta sequence for improved spectral sampling. |
| IIa-prime | `12 x 12 x 4` | `18 x 18 x 6` | `3.375x`; `Delta k_max = 0.103625 Angstrom^-1` | Applies the same 1.5-fold directional refinement while retaining fine sampling along the longer c repeat. |
| IIb | `12 x 12 x 4` | `18 x 18 x 6` | `3.375x`; `Delta k_max = 0.103794 Angstrom^-1` | Applies the same 1.5-fold directional refinement while retaining fine sampling along the longer c repeat. |

The DOS/PDOS meshes are deliberately denser than the authoritative Stage 03 total-energy meshes to improve spectral sampling quality. Their use does not redefine the converged polymorph energies: energy convergence was established by Stage 02 and the authoritative five-phase energies remain those of Stage 03.

## Band-Structure KPOINTS

Conventional band structures require ordered trajectories between labelled high-symmetry points. The existing Stage 04 band files therefore use VASP reciprocal-coordinate Line-mode rather than automatic three-dimensional Gamma meshes. Inspection confirms the following implemented paths and segment densities:

| Phase | Space Group | SeeK-path Type | Canonical Path | Points per Segment |
| --- | --- | --- | --- | ---: |
| Spinel | `Fd-3m` (No. 227) | `cF2` | `GAMMA-X-U | K-GAMMA-L-W-X` | 40 |
| Alpha1 | `R3m` (No. 160) | `hR1` | `GAMMA-T-H_2 | H_0-L-GAMMA-S_0 | S_2-F-GAMMA` | 40 |
| Beta | `P3m1` (No. 156) | `hP2` | `GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K` | 40 |
| IIa-prime | `P-3m1` (No. 164) | `hP2` | `GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K` | 40 |
| IIb | `P6_3mc` (No. 186) | `hP2` | `GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K` | 40 |

## SeeK-path Validation

The band paths follow the HPKOT/SeeK-path convention and were independently validated using the official SeeK-path package version 2.2.1 with spglib 2.7.0 and the default `symprec = 1e-5`.

SeeK-path internally constructs a standard crystallographic representation. For this validation, each canonical endpoint was first interpreted as a Cartesian reciprocal-space vector using the SeeK-path reciprocal primitive lattice and rotation, then expressed in the reciprocal basis of the original, unchanged frozen POSCAR. The maximum Cartesian endpoint mismatch across all five existing KPOINTS files was `1.386 x 10^-10 Angstrom^-1`. This is numerically negligible and confirms both the path definitions and the basis transformations.

SeeK-path did not directly write these files. No standardised POSCAR was substituted into the VASP inputs: only the canonical path coordinates were transformed back to the original POSCAR basis.

## POSCAR / CHGCAR Compatibility

Changing KPOINTS does not by itself invalidate a CHGCAR when the real-space cell and structure used to interpret that charge density are preserved. Changing the lattice vectors, cell representation, or associated real-space grid could invalidate that compatibility even if the replacement cell were crystallographically equivalent.

The Stage 04 design therefore keeps each band and DOS/PDOS POSCAR identical to its corresponding frozen Stage 03 POSCAR. Direct inspection also confirms that each current Stage 04 CHGCAR is byte-identical to the corresponding Stage 03 CHGCAR. The band endpoints are expressed in that original POSCAR's reciprocal basis, preserving the exact structure/charge-density lineage required for the non-self-consistent `ICHARG = 11` calculations. This does not imply that CHGCAR files may be transferred arbitrarily between differently represented cells.

## Symmetry-Tolerance Note

At unusually strict symmetry tolerances below approximately `2 x 10^-6`, tiny relaxed-coordinate distortions cause alpha1 and IIb to lose their nominal higher symmetry in automated detection. Both are correctly identified as `R3m` and `P6_3mc`, respectively, at the official SeeK-path default `symprec = 1e-5`. This numerical sensitivity does not affect the validated Stage 04 band-path selection.

## Final Selection

The accepted Stage 04 strategy is:

- **DOS/PDOS:** dense, uniform, unshifted Gamma-centred meshes: `6 x 6 x 6` for spinel and `18 x 18 x 6` for alpha1, beta, IIa-prime, and IIb.
- **Band structures:** official HPKOT/SeeK-path high-symmetry Line-mode paths, sampled with 40 points per segment and represented in the reciprocal basis of each original frozen POSCAR.

The existing inputs meet these selections. No further KPOINTS validation is required before running Stage 04.
